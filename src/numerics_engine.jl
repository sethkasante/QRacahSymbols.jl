# src/Numerics.jl

# const BIG_PI = big(π)

struct NumericSU2kModel{T <: AbstractFloat}
    k::Int
    logqnfact::Vector{T}
end

function NumericSU2kModel(k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
    tab::Vector{T} = get!(LOGQFACT_CACHE, (k, T, prec)) do
        logqnfact_table(k, T, prec) 
    end
    return NumericSU2kModel{T}(k, tab)
end

function qinteger_num(n::Int, k::Int, T::Type)
    θ = big(π) / T(k+2)
    return sin(n * θ) / sin(θ)
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
            θ = big(π) / BigFloat(N)
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

@inline function log_qΔ(j1, j2, j3, tab::Vector{T})::T where {T}
    a = Int(j1 + j2 - j3); b = Int(j1 - j2 + j3); c = Int(-j1 + j2 + j3); d = Int(j1 + j2 + j3)
    return (tab[a+1] + tab[b+1] + tab[c+1] - tab[d+2]) / 2 
end

@inline function logqtri_coeffs(j1, j2, j3, j4, j5, j6, tab::Vector{T})::T where {T}
    return log_qΔ(j1, j2, j3, tab) + log_qΔ(j1, j5, j6, tab) + log_qΔ(j2, j4, j6, tab) + log_qΔ(j3, j4, j5, tab)
end

@inline function log_racah_summand(z, α1, α2, α3, α4, β1, β2, β3, tab::Vector{T})::T where {T}
    lognum = tab[z+2]
    logden = tab[z-α1+1] + tab[z-α2+1] + tab[z-α3+1] + tab[z-α4+1] + tab[β1-z+1] + tab[β2-z+1] + tab[β3-z+1]
    return lognum - logden
end

function _qracah6j_stable(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    # if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
    #     return zero(T)
    # end
    #admissible checks in main file
    
    table = model.logqnfact
    logT = logqtri_coeffs(j1, j2, j3, j4, j5, j6, table)

    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6); α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)

    zrange = max(α1, α2, α3, α4):min(β1, β2, β3, model.k) 
    
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

# Disambiguated Wrappers
qracah6j_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin) where {T} =
    _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)

# function qracah6j_numeric(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
#     model = NumericSU2kModel(k; T=T, prec=prec)
#     return _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
# end
function qracah6j_numeric(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
    if T === BigFloat
        return setprecision(BigFloat, prec) do
            model = NumericSU2kModel(k; T=T, prec=prec)
            _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
        end
    else
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    end
end



# ============================================================
# High-Performance Numeric 3j Evaluator
# ============================================================

@inline function log_q3j_prefactor(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, tab::Vector{T})::T where {T}
    log_delta = log_qΔ(j1, j2, j3, tab)
    log_facts = tab[Int(j1+m1)+1] + tab[Int(j1-m1)+1] + 
                tab[Int(j2+m2)+1] + tab[Int(j2-m2)+1] + 
                tab[Int(j3-m1-m2)+1] + tab[Int(j3+m1+m2)+1]
                
    # log(sqrt(Delta^2 * facts)) = log_delta + 0.5 * log_facts
    return log_delta + 0.5 * log_facts
end

@inline function log_q3j_summand(z::Int, α1::Int, α2::Int, β1::Int, β2::Int, β3::Int, tab::Vector{T})::T where {T}
    # Numerator is 1, so log(1) = 0. We just negate the sum of the denominator logs.
    return -(tab[z+1] + tab[α1+z+1] + tab[α2+z+1] + tab[β1-z+1] + tab[β2-z+1] + tab[β3-z+1])
end

function _qracah3j_stable(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin)::T where {T}
    #admissible checks performed in the main file
    table = model.logqnfact

    logT = log_q3j_prefactor(j1, j2, j3, m1, m2, table)

    α1 = Int(j3 - j2 + m1) 
    α2 = Int(j3 - j1 - m2)
    β1 = Int(j1 + j2 - j3)
    β2 = Int(j1 - m1)
    β3 = Int(j2 + m2)

    zrange = max(-α1, -α2, 0):min(β1, β2, β3, model.k) 
    
    logmax = typemin(T)
    @inbounds for z in zrange
        logTsz = logT + log_q3j_summand(z, α1, α2, β1, β2, β3, table)
        logmax = max(logmax, logTsz)
    end

    res_scaled = zero(T)
    phase_offset = α1 - α2
    
    @inbounds for z in zrange
        logTsz = logT + log_q3j_summand(z, α1, α2, β1, β2, β3, table)
        val = exp(logTsz - logmax)
        res_scaled += isodd(z + phase_offset) ? -val : val
    end

    return exp(logmax) * res_scaled
end

qracah3j_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin) where {T} =
    _qracah3j_stable(model, j1, j2, j3, m1, m2)

# function qracah3j_numeric(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
#     model = NumericSU2kModel(k; T=T, prec=prec)
#     return _qracah3j_stable(model, j1, j2, j3, m1, m2)
# end

function qracah3j_numeric(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, k::Int; T::Type{<:AbstractFloat}=Float64, prec::Int=256)
    if T === BigFloat
        return setprecision(BigFloat, prec) do
            model = NumericSU2kModel(k; T=T, prec=prec)
            _qracah3j_stable(model, j1, j2, j3, m1, m2)
        end
    else
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _qracah3j_stable(model, j1, j2, j3, m1, m2)
    end
end