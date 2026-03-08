# ============================================================
# 1. Quantum Integers [n] and Dimensions [2j+1]_q
# ============================================================

"""
    qint(n::Int; mode=:generic)
    qint(n::Int, k::Int; mode=:exact, T=Float64)

Returns a representation of the quantum integer [n]_q. 
The `mode` determines the choice of representation:
- `:generic` : Returns a symbolic `CycloMonomial` factorization (requires only `n`).
- `:exact`   : Returns the exact algebraic representation in the cyclotomic field Q(ζ).
- `:numeric` : Returns the numerical floating-point value using the stable sine formula.
"""
function qint end

# 1. Symbolic Mode (Requires only n)
function qint(n::Int; mode=:generic)
    @assert mode == :generic "Numeric and exact modes require level k. Use qint(n, k; mode=...)"
    n == 0 && return CycloMonomial(0, 0, Int[]) 
    n == 1 && return CycloMonomial(1, 0, Int[]) 
    
    exps = zeros(Int, n)
    for d in 2:n
        if n % d == 0
            exps[d] = 1
        end
    end
    return CycloMonomial(1, 0, exps)
end

# 2. Exact and Numeric Modes (Requires n and k)
function qint(n::Int, k::Int; mode=:exact, T::Type{<:AbstractFloat}=Float64)
    if mode == :exact
        n == 0 && return Nemo.QQFieldElem(0) # Generic zero fallback
        n == 1 && return Nemo.QQFieldElem(1)
        
        # Smart routing: Use cache if it exists, otherwise compute on the fly
        if haskey(EXACT_MODEL_CACHE, k)
            model = EXACT_MODEL_CACHE[k]
            return model.q_facts[n+1] * inv(model.q_facts[n])
        else
            # Lightweight on-the-fly calculation
            K, z = Nemo.cyclotomic_field(2 * (k + 2), "ζ")
            n == 0 && return K(0)
            n == 1 && return K(1)
            
            # [n]_q = (z^n - z^-n) / (z - z^-1)
            z_inv = inv(z)
            return (z^n - z_inv^n) * inv(z - z_inv)
        end
        
    elseif mode == :numeric
        n == 0 && return zero(T)
        n == 1 && return one(T)
        θ = T(π) / T(k + 2)
        return sin(T(n) * θ) / sin(θ)
        
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end


function qdim_symb(j::Spin)
    n = Int(2j + 1)
    buf = SymbolicBuffer(n + 1)
    add_qint!(buf, n, 1)
    return snapshot(buf)
end

function qdim_exact(model::ExactSU2kModel, j::Spin)
    n = Int(2j + 1)
    return model.q_facts[n+1] * inv(model.q_facts[n])
end

function qdim_numeric(j::Spin, model::NumericSU2kModel{T})::T where {T}
    n = Int(2j + 1)
    return exp(model.logqnfact[n+1] - model.logqnfact[n])
end

# ============================================================
# 2. R-Matrix (Braiding)
# ============================================================

function rmatrix_symb(j1::Spin, j2::Spin, j3::Spin)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    return CycloMonomial(s, phase_exp, Int[])
end

function rmatrix_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    res = model.z^phase_exp
    return s == 1 ? res : -res
end

function rmatrix_numeric(j1::Spin, j2::Spin, j3::Spin, k::Int; T::Type{<:AbstractFloat}=Float64)
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    return s * cispi(T(phase_exp) / T(k + 2))
end

# ============================================================
# 3. F-Symbol (Fusion / Normalized 6j)
# ============================================================

function fsymbol_generic(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq.sign == 0 && return res_6j
    
    # 1. Update the prefactor using the SymbolicBuffer to avoid CycloMonomial multiplications
    cap = length(res_6j.pref_sq.exps) + Int(2 * max(j3, j6) + 2)
    buf = SymbolicBuffer(cap)
    
    # Initialize buffer with the existing 6j prefactor
    buf.sign = res_6j.pref_sq.sign
    buf.z_pow = res_6j.pref_sq.z_pow
    copyto!(buf.exps, 1, res_6j.pref_sq.exps, 1, length(res_6j.pref_sq.exps))
    
    # Add dimensions [2j3+1]_q and [2j6+1]_q
    add_qint!(buf, Int(2j3 + 1), 1)
    add_qint!(buf, Int(2j6 + 1), 1)
    
    new_pref_sq = snapshot(buf)
    
    # 2. Apply overall phase shift to the Racah sum series
    phase_s = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    if phase_s == -1
        new_series = [CycloMonomial(-term.sign, term.z_pow, copy(term.exps)) for term in res_6j.series]
        return GenericResult(new_pref_sq, new_series)
    end
    
    return GenericResult(new_pref_sq, res_6j.series)
end

function fsymbol_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq == 0 && return res_6j
    
    new_pref_sq = res_6j.pref_sq * qdim_exact(model, j3) * qdim_exact(model, j6)
    phase_s = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    new_sum = phase_s == 1 ? res_6j.sum_cf : -res_6j.sum_cf
    
    return ExactResult(model.k, new_pref_sq, new_sum)
end

function fsymbol_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    phase = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    return phase * sqrt(qdim_numeric(j3, model) * qdim_numeric(j6, model)) * val_6j
end

# ============================================================
# 4. G-Symbol (Tetrahedral Weight)
# ============================================================

function gsymbol_generic(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_generic(j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq.sign == 0 && return res_6j
    
    cap = length(res_6j.pref_sq.exps) + Int(2 * max(j1, j2, j3, j4, j5, j6) + 2)
    buf = SymbolicBuffer(cap)
    
    buf.sign = res_6j.pref_sq.sign
    buf.z_pow = res_6j.pref_sq.z_pow
    copyto!(buf.exps, 1, res_6j.pref_sq.exps, 1, length(res_6j.pref_sq.exps))
    
    add_qint!(buf, Int(2j1 + 1), 1); add_qint!(buf, Int(2j2 + 1), 1); add_qint!(buf, Int(2j3 + 1), 1)
    add_qint!(buf, Int(2j4 + 1), 1); add_qint!(buf, Int(2j5 + 1), 1); add_qint!(buf, Int(2j6 + 1), 1)
    
    return GenericResult(snapshot(buf), res_6j.series)
end

function gsymbol_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
    res_6j.pref_sq == 0 && return res_6j
    
    dim_prod = qdim_exact(model, j1) * qdim_exact(model, j2) * qdim_exact(model, j3) * qdim_exact(model, j4) * qdim_exact(model, j5) * qdim_exact(model, j6)
    return ExactResult(model.k, res_6j.pref_sq * dim_prod, res_6j.sum_cf)
end

function gsymbol_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    val_6j = _qracah6j_stable(model, j1, j2, j3, j4, j5, j6)
    dims_prod = qdim_numeric(j1, model) * qdim_numeric(j2, model) * qdim_numeric(j3, model) * qdim_numeric(j4, model) * qdim_numeric(j5, model) * qdim_numeric(j6, model)
    return val_6j * sqrt(dims_prod)
end