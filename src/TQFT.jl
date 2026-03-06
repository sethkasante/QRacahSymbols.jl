# src/tqft.jl

# ============================================================
# 1. Quantum Dimensions [2j+1]_q
# ============================================================
function qdim_symb(j::Spin)
    n = Int(2j + 1)
    return qfactorial_symb(n) / qfactorial_symb(n - 1)
end

function qdim_exact(model::ExactSU2kModel, j::Spin)
    return evaluate_cyclofield(qdim_symb(j), model)
end

function qdim_numeric(j::Spin, model::NumericSU2kModel{T})::T where {T}
    n = Int(2j + 1)
    return exp(model.logqnfact[n+1] - model.logqnfact[n])
end

# ============================================================
# 2. R-Matrix (Braiding)
# ============================================================
function rmatrix_symb(j1::Spin, j2::Spin, j3::Spin)
    !δ(j1, j2, j3) && return CycloMonomial(0, 0, Int[])
    phase_exp = Int(j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(Int(j1 + j2 - j3)) ? 1 : -1
    return CycloMonomial(s, phase_exp, Int[])
end

function rmatrix_exact(model::ExactSU2kModel, j1::Spin, j2::Spin, j3::Spin)
    return evaluate_cyclofield(rmatrix_symb(j1, j2, j3), model)
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
    
    new_pref_sq = res_6j.pref_sq * qdim_symb(j3) * qdim_symb(j6)
    phase_s = iseven(Int(j1 + j2 + j4 + j5)) ? 1 : -1
    
    new_series = phase_s == -1 ? [CycloMonomial(-term.sign, term.z_pow, term.exps) for term in res_6j.series] : res_6j.series
    return GenericResult(new_pref_sq, new_series)
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
    
    dim_prod = qdim_symb(j1) * qdim_symb(j2) * qdim_symb(j3) * qdim_symb(j4) * qdim_symb(j5) * qdim_symb(j6)
    return GenericResult(res_6j.pref_sq * dim_prod, res_6j.series)
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