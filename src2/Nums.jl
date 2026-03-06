# src/Numerics.jl


#Major #TODO: improve type stability.  

# Precomputing BigFloat Pi to avoid repeated allocations in the BigFloat path
const BIG_PI = big(π)

# Cache key: (k, Type, precision) -> Vector. 
# Storing as `Any` in the LRU but type-asserting on retrieval ensures type-stability.
const LOGQFACT_CACHE = LRU{Tuple{Int, DataType, Int}, Any}(maxsize = 10240)


#TODO: Is it better to embed NumericSU2kModel function inside struct? 
# Parameterize the model to accept Float64 OR BigFloat
struct NumericSU2kModel{T <: AbstractFloat}
    k::Int
    logqnfact::Vector{T}
end

# Default to Float64 for blistering speed, opt-in to BigFloat for extreme precision
function NumericSU2kModel(k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
    # Type assertion ::Vector{T} guarantees the compiler knows exactly what comes out
    tab::Vector{T} = get!(LOGQFACT_CACHE, (k, T, prec)) do
        logqnfact_table(k, T, prec) 
    end
    return NumericSU2kModel{T}(k, tab)
end

function logqn_table(k::Int, T::Type, prec::Int)::Vector{T}
    N = k + 2
    
    if T === Float64
        θ = pi / N
        logsinθ = log(sin(θ))
        tab = Vector{T}(undef, N)
        tab[1] = zero(T)
        half = N ÷ 2
        @inbounds for n in 1:half
            val = log(sin(n * θ)) - logsinθ
            tab[n+1] = val
            tab[N + 1 - n] = val
        end
        return tab
    else
        setprecision(BigFloat, prec) do
            θ = BIG_PI / BigFloat(N)
            logsinθ = log(sin(θ))
            tab = Vector{T}(undef, N)
            tab[1] = zero(T)
            half = N ÷ 2
            @inbounds for n in 1:half
                val = log(sin(n * θ)) - logsinθ
                tab[n+1] = val
                tab[N + 1 - n] = val
            end
            return tab
        end
    end
end

function logqnfact_table(k::Int, T::Type, prec::Int)::Vector{T}
    logqn = logqn_table(k, T, prec)
    tab = Vector{T}(undef, k + 2) 
    tab[1] = logqn[1]  
    @inbounds for n in 2:k+2
        tab[n] = tab[n-1] + logqn[n]
    end
    return tab
end

# ============================================================
# Core Log Building Blocks (Strictly Generic)
# ============================================================

@inline function log_qΔ(j1, j2, j3, tab::Vector{T})::T where {T}
    a = Int(j1 + j2 - j3)
    b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3)
    d = Int(j1 + j2 + j3)
    # Division by an integer literal dynamically preserves T
    return (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+2]) / 2 
end

@inline function logqtri_coeffs(j1, j2, j3, j4, j5, j6, tab::Vector{T})::T where {T}
    return log_qΔ(j1, j2, j3, tab) + 
           log_qΔ(j1, j5, j6, tab) +
           log_qΔ(j2, j4, j6, tab) + 
           log_qΔ(j3, j4, j5, tab)
end

@inline function log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, tab::Vector{T})::T where {T}
    lognum = tab[z+2]
    logden = tab[z-α1+1] + tab[z-α2+1] + tab[z-α3+1] +
             tab[z-α4+1] + tab[β1-z+1] + tab[β2-z+1] + tab[β3-z+1]
    return lognum - logden
end

# ============================================================
# High-Performance Numeric 6j Evaluator
# ============================================================

"""
    _qracah6j_stable(model::NumericSU2kModel{T}, j1, j2, j3, j4, j5, j6) where {T}

Computes the quantum 6j symbol. Return type scales with the Model parameter T.
"""
function _qracah6j_stable(model::NumericSU2kModel{T}, j1, j2, j3, j4, j5, j6)::T where {T}
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return zero(T)
    end
    
    table = model.logqnfact

    logT = logqtri_coeffs(j1, j2, j3, j4, j5, j6, table)

    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)

    zrange = max(α1, α2, α3, α4):min(β1, β2, β3, model.k) 
    
    # Use typemin(T) instead of -Inf to guarantee type stability
    logmax = typemin(T)
    @inbounds for z in zrange
        logTsz = logT + log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, table)
        logmax = max(logmax, logTsz)
    end

    res_scaled = zero(T)
    @inbounds for z in zrange
        logTsz = logT + log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, table)
        val = exp(logTsz - logmax)
        res_scaled += iseven(z) ? val : -val
    end

    return exp(logmax) * res_scaled
end

# Model evaluation wrappers
qracah6j_numeric(model::NumericSU2kModel{T}, j1, j2, j3, j4, j5, j6) where {T} =
    _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)

"""
    qracah6j_numeric(j1, j2, j3, j4, j5, j6, k::Int; T::Type=Float64, prec::Int=256)

Standalone wrapper for evaluating the numeric 6j symbol for a specific `k`.
Creates or retrieves a NumericSU2kModel from the cache automatically.
"""
function qracah6j_numeric(k::Int, j1, j2, j3, j4, j5, j6; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
    model = NumericSU2kModel(k; T=T, prec=prec)
    return _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
end


#TODO: Fix types of inputs?: sometimes qracah6j_numeric can be error with ambiguous inputs