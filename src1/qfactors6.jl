module SymbolicQRacah

using Nemo
using LRUCache

export CycloMonomial, ExactSU2kModel, qfactorial_symb, q6jseries_symb
export exact_qracah6j, qracah6j_symb, qracah6jsq_symb

# ============================================================
# 1. Symbolic Data Structures (Quantum Primes)
# ============================================================

"""
    CycloMonomial

Represents exactly a product of cyclotomic polynomials evaluated at q:
Value = sign * z^(z_pow) * Π (Φ_d(z^2))^(exps[d])

Here, z = q^{1/2} = e^{i π / (k+2)}.
"""
struct CycloMonomial
    sign::Int
    z_pow::Int           # Handles the integer powers of z (which are fractional powers of q)
    exps::Dict{Int, Int} # Handles Π Φ_d(z^2)^{exps[d]}
end

import Base: *, /

function *(a::CycloMonomial, b::CycloMonomial)
    small, large = length(a.exps) ≤ length(b.exps) ? (a,b) : (b,a)
    exps = copy(large.exps)
    for (d, e) in small.exps
        val = get(exps, d, 0) + e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    return CycloMonomial(a.sign * b.sign, a.z_pow + b.z_pow, exps)
end

function /(a::CycloMonomial, b::CycloMonomial)
    exps = copy(a.exps)
    for (d, e) in b.exps
        val = get(exps, d, 0) - e
        if val == 0
            delete!(exps, d)
        else
            exps[d] = val
        end
    end
    return CycloMonomial(a.sign * b.sign, a.z_pow - b.z_pow, exps)
end

# ============================================================
# Human-Readable Output (Pretty Printing)
# ============================================================

# Unicode mappings for clean subscripts and superscripts
const SUBSCRIPTS = Dict('0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', 
                        '5'=>'₅', '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉')
const SUPERSCRIPTS = Dict('0'=>'⁰', '1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', 
                          '5'=>'⁵', '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '-'=>'⁻')

function to_subscript(n::Int)
    return map(c -> SUBSCRIPTS[c], string(n))
end

function to_superscript(n::Int)
    return map(c -> SUPERSCRIPTS[c], string(n))
end

# Overload Base.show to format the REPL output
function Base.show(io::IO, M::CycloMonomial)
    # Edge case for exact zero
    if M.sign == 0
        print(io, "0")
        return
    end

    parts = String[]
    
    # 1. z power (z = q^{1/2})
    if M.z_pow != 0
        if M.z_pow == 1
            push!(parts, "z")
        else
            push!(parts, "z" * to_superscript(M.z_pow))
        end
    end
    
    # 2. Cyclotomic polynomials Φ_d
    # Sort keys to ensure deterministic, mathematically neat ordering
    for d in sort(collect(keys(M.exps)))
        e = M.exps[d]
        e == 0 && continue
        
        base_str = "Φ" * to_subscript(d)
        if e == 1
            push!(parts, base_str)
        else
            push!(parts, base_str * to_superscript(e))
        end
    end
    
    # 3. Assembly
    if isempty(parts)
        # If there are no z powers or Φ terms, it's just the constant 1 or -1
        print(io, M.sign == -1 ? "-1" : "1")
    else
        # Prefix the sign, and join the terms with a space
        prefix = M.sign == -1 ? "-" : ""
        print(io, prefix * join(parts, " "))
    end
end


# ============================================================
# 2. Exact Algebraic Model (Nemo.jl)
# ============================================================

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     # Primitive 2N-th root of unity (z = q^{1/2})
    Phi_cache::Dict{Int, nf_elem}  # Cached evaluations of Φ_d(z^2)
end

function ExactSU2kModel(k::Int)
    N = k + 2
    # z is a primitive 2N-th root of unity: z = exp(i * pi / N)
    K, z = cyclotomic_field(2N, "ζ") 
    
    Zx, x = polynomial_ring(ZZ, "x") 
    
    # We cache evaluations of Φ_d at z^2 (which represents q)
    # The maximum value inside a factorial is 2k (which is < 2N)
    max_d = 2 * N 
    Phi_cache = Dict{Int, nf_elem}()
    
    z_sq = z^2
    for d in 2:max_d
        poly = cyclotomic(d, x)
        Phi_cache[d] = evaluate(poly, z_sq)  
    end

    return ExactSU2kModel(k, K, z, Phi_cache)
end

# ============================================================
# 4. Symbolic Math Engine
# ============================================================

function qfactorial_symb(n::Int)
    if n == 0
        return CycloMonomial(1, 0, Dict{Int, Int}())
    end

    exps = Dict{Int, Int}()
    for d in 2:n
        power = div(n, d) 
        if power > 0
            exps[d] = power
        end
    end
    
    # Power of z = q^{1/2}. q^{-n(n-1)/4} = z^{-n(n-1)/2}
    z_pow = - (n * (n - 1)) ÷ 2
    
    return CycloMonomial(1, z_pow, exps)
end

function qdelta2_symb(a, b, c)
    num = qfactorial_symb(Int(a+b-c)) * qfactorial_symb(Int(a-b+c)) * qfactorial_symb(Int(-a+b+c))
    den = qfactorial_symb(Int(a+b+c+1))
    return num / den
end

function qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    return qdelta2_symb(j1, j2, j3) * qdelta2_symb(j1, j5, j6) * qdelta2_symb(j2, j4, j6) * qdelta2_symb(j3, j4, j5)
end

function q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3)
    num = qfactorial_symb(z+1)
    den = qfactorial_symb(z-α1) * qfactorial_symb(z-α2) * qfactorial_symb(z-α3) *
          qfactorial_symb(z-α4) * qfactorial_symb(β1-z) * qfactorial_symb(β2-z) *
          qfactorial_symb(β3-z) 
    
    res = num / den
    if isodd(z)
        return CycloMonomial(-res.sign, res.z_pow, res.exps)
    else
        return res
    end
end

function q6jseries_symb(j1, j2, j3, j4, j5, j6)::Vector{CycloMonomial}
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    S_z = CycloMonomial[]
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3)
    @inbounds for z in zrange
        push!(S_z, q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3)) 
    end
    return S_z
end


# ============================================================
# 5. Field Evaluation & BigFloat Bridge
# ============================================================

"""
Maps a purely symbolic CycloMonomial exactly into the Nemo number field.
"""
function evaluate_nf(M::CycloMonomial, model::ExactSU2kModel)
    K = model.K
    res = K(M.sign) * model.z^M.z_pow
    for (d, e) in M.exps
        res *= (model.Phi_cache[d])^e
    end
    return res
end


# ============================================================
# 6. Final Evaluator API
# ============================================================

@inline ishalfInt(j)::Bool = isinteger(2*j) && j ≥ 0
@inline δ(j1, j2, j3)::Bool = isinteger(j1+j2+j3) && (abs(j1-j2) <= j3 <= j1+j2)
@inline qδ(j1, j2, j3, k)::Bool = δ(j1, j2, j3) && (j1 + j2 + j3) <= k 

@inline function qδtet(j1, j2, j3, j4, j5, j6, k)
    return ishalfInt(j1) && ishalfInt(j2) && ishalfInt(j3) && 
           ishalfInt(j4) && ishalfInt(j5) && ishalfInt(j6) && 
           qδ(j1, j2, j3, k) && qδ(j1, j5, j6, k) && 
           qδ(j2, j4, j6, k) && qδ(j3, j4, j5, k)
end

#exact symbolic evaluation 
#for now we keep the exact evaluation in two terms first is sqrt of triangle coefficient and second is the sum. Is there a better way to present it?

function qracah6j_symb(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return zero(BigFloat)
    end
    
    # 1. Exact Summation in the Field
    S_z = q6jseries_symb(j1, j2, j3, j4, j5, j6)
    sum_nf = model.K(0)
    for M in S_z
        sum_nf += evaluate_nf(M, model)
    end
    
    # 2. Exact Prefactor ^ 2 in the Field
    Pref2_symb = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    Pref2_nf = evaluate_nf(Pref2_symb, model)
    return Pref2_nf, sum_nf
end

function exact_qracah6j(j1, j2, j3, j4, j5, j6,k::Int)
    model = ExactSU2kModel(k)
    qracah6j_symb(model, j1, j2, j3, j4, j5, j6)
end

# The square of the racah symbol. How do we keep track of the sign? In case we want to take sqrt of this? 
function qracah6jsq_symb(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    Pref2_nf, sum_nf = qracah6j_symb(model, j1, j2, j3, j4, j5, j6)
    return Pref2_nf * sum_nf^2
    #What about the sign? 
end


# for when numeric is needed : not necessary here but we keep it anyways. 
"""
Safely evaluates an exact number field element to a Complex{BigFloat} using Horner's method.
"""
function evaluate_horner(ev::nf_elem, root_val)
    d = degree(parent(ev)) - 1
    d < 0 && return zero(Complex{BigFloat})
    
    res = BigFloat(coeff(ev, d))
    for i in (d-1):-2:0
        c = BigFloat(coeff(ev, i))
        res = res * root_val + c
    end
    return res
end

end # module



























# src/ExactAlgebra.jl

using Nemo

# ============================================================
# Exact Algebraic Model (Nemo.jl)
# ============================================================

struct ExactSU2kModel
    k::Int
    K::AnticNumberField            
    z::nf_elem                     
    q_facts::Vector{nf_elem}       # Precomputed exact [n]_q! in the field
end

"""
    ExactSU2kModel(k::Int)

Instantiates the Cyclotomic field for level k and precomputes the 
exact quantum factorials natively in the field.
"""
function ExactSU2kModel(k::Int)
    N = k + 2
    K, z = cyclotomic_field(2N, "ζ") # z = q^1/2 = exp(i π/(k+2)) 
    
    # Precompute exact quantum factorials [n]_q! natively
    # Array size is k+3 to safely store indices for [0]_q! through [k+2]_q!
    q_facts = Vector{nf_elem}(undef, k + 3)
    q_facts[1] = K(1) # [0]_q! = 1
    
    # Efficiently build [n]_q = (z^n - z^-n) / (z - z^-1)
    z_inv = inv(z)
    den_inv = inv(z - z_inv)
    z_n = z
    z_inv_n = z_inv
    
    @inbounds for n in 1:(k+2)
        q_int = (z_n - z_inv_n) * den_inv
        q_facts[n+1] = q_facts[n] * q_int
        
        z_n *= z
        z_inv_n *= z_inv
    end

    return ExactSU2kModel(k, K, z, q_facts)
end

# ============================================================
# High-Speed Exact Algebraic 6j Evaluator
# ============================================================

@inline function q_delta2_exact(model::ExactSU2kModel, j1, j2, j3)
    a = Int(j1 + j2 - j3)
    b = Int(j1 - j2 + j3)
    c = Int(-j1 + j2 + j3)
    d = Int(j1 + j2 + j3)
    
    # Instant memory lookup using precomputed field factorials
    num = model.q_facts[a+1] * model.q_facts[b+1] * model.q_facts[c+1]
    den = model.q_facts[d+2]
    return num * inv(den)
end

@inline function qtricoeff2_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    return q_delta2_exact(model, j1, j2, j3) * q_delta2_exact(model, j1, j5, j6) * q_delta2_exact(model, j2, j4, j6) * q_delta2_exact(model, j3, j4, j5)
end

function q6jseries_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3, model.k) 
    sum_cf = model.K(0)

    @inbounds for z in zrange
        num = model.q_facts[z+2]
        den = model.q_facts[z - α1 + 1] * model.q_facts[z - α2 + 1] * model.q_facts[z - α3 + 1] * model.q_facts[z - α4 + 1] * model.q_facts[β1 - z + 1] * model.q_facts[β2 - z + 1] * model.q_facts[β3 - z + 1]
              
        term = num * inv(den)
        sum_cf += iseven(z) ? term : -term
    end
    
    return sum_cf
end

function qracah6j_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return ExactResult(model.k, model.K(0), model.K(0))
    end
    
    Tc2 = qtricoeff2_exact(model, j1, j2, j3, j4, j5, j6)
    Sum_cf = q6jseries_exact(model, j1, j2, j3, j4, j5, j6)
    
    return ExactResult(model.k, Tc2, Sum_cf)
end

# ============================================================
# Float Projection (Horner's Method)
# ============================================================

function evaluate_exact(res::ExactResult, prec::Int=256)
    if res.pref_sq == 0 && res.sum_cf == 0
        return zero(BigFloat)
    end
    
    setprecision(BigFloat, prec) do
        N = res.k + 2
        target_z = cispi(big"1.0" / N)
        
        sum_bf = horner_eval(res.sum_cf, target_z)
        pref2_bf = horner_eval(res.pref_sq, target_z)
        
        return real(sqrt(pref2_bf) * sum_bf)
    end
end

function horner_eval(ev::nf_elem, root_val::Complex{BigFloat})
    d = degree(parent(ev)) - 1
    d < 0 && return zero(Complex{BigFloat})
    
    res = Complex{BigFloat}(coeff(ev, d))
    for i in (d-1):-1:0
        c = coeff(ev, i)
        val = BigFloat(Nemo.numerator(c)) / BigFloat(Nemo.denominator(c))
        res = res * root_val + val
    end
    return res
end

























# src/ExactAlgebra.jl


# Global Polynomial Cache (Independent of k)
# ============================================================

# Create the generic integer polynomial ring exactly once globally
const ZX_RING, ZX_VAR = polynomial_ring(ZZ, "x")

# Cache the abstract cyclotomic polynomials Φ_d(x)
const CYCLO_POLY_CACHE = LRU{Int, elem_type(ZX_RING)}(maxsize = 1024)


"""
    get_cyclo_poly(d::Int)

Retrieves the abstract integer cyclotomic polynomial Φ_d(x) from the global cache,
or computes and caches it if it doesn't exist.
"""
@inline function get_cyclo_poly(d::Int)
    return get!(CYCLO_POLY_CACHE, d) do
        cyclotomic(d, ZX_VAR)
    end
end


# ============================================================
# Exact Algebraic Model using (Nemo.jl)
# ============================================================

function ExactSU2kModel(k::Int)
    N = k + 2
    # z is a primitive 2N-th root of unity: z = exp(i * pi / N)
    K, z = cyclotomic_field(2N, "ζ") 
    
    # We cache evaluations of Φ_d at z^2 (which represents q)
    max_d = 2 * N 
    
    # Pre-allocate a dense vector for O(1) field evaluations during Racah sums
    Phi_eval = Vector{nf_elem}(undef, max_d)
    
    z_sq = z^2
    @inbounds for d in 1:max_d
        if d == 1
            Phi_eval[d] = z_sq - K(1) # Φ_1(q) = q - 1
        else
            # 1. Grab the abstract polynomial from the GLOBAL cache
            poly = get_cyclo_poly(d)
            # 2. Evaluate it in the LOCAL number field
            Phi_eval[d] = evaluate(poly, z_sq)  
        end
    end

    # Return the model with its local field evaluation cache
    return ExactSU2kModel(k, K, z, Phi_eval)
end


function ExactSU2kModel(k::Int)
    # ... existing K, z setup ...
    
    # 1. Precompute exact quantum integers [n]_q in the field
    # [n]_q = (z^(n) - 1) / (z - 1)
    q_ints = Vector{nf_elem}(undef, k + 2)
    q_ints[1] = K(0)
    if k >= 0 
        q_ints[2] = K(1) 
    end
    
    # Efficiently build field integers
    z_sq = z^2
    z_2n = z_sq^2
    denom_inv = inv(z_sq - K(1))
    for n in 3:(k+2)
        q_ints[n] = (z_2n - K(1)) * denom_inv
        z_2n *= z_sq
    end

    # 2. Precompute exact factorials [n]_q!
    q_facts = Vector{nf_elem}(undef, k + 2)
    q_facts[1] = K(1) # [0]! = 1
    for n in 2:(k+2)
        q_facts[n] = q_facts[n-1] * q_ints[n]
    end

    return ExactSU2kModel(k, K, z, q_facts) # Throw away Phi_cache entirely!
end

# ============================================================
# Field Evaluation in Cyclotomic Fields 
# ============================================================

"""
Maps a purely symbolic CycloMonomial exactly into the Cyclotomic number field.
"""
function evaluate_cyclofield(M::CycloMonomial, model::ExactSU2kModel)
    K = model.K
    res = K(M.sign) 

    # Safe exponentiation for field generators
    if M.z_pow != 0
        base_z = M.z_pow > 0 ? model.z : inv(model.z)
        res *= base_z^abs(M.z_pow)
    end
    
    # Iterate over the Vector indices (which represent the cyclotomic base d)
    @inbounds for d in 1:length(M.exps)
        e = M.exps[d]
        if e != 0
            term = model.Phi_eval[d]
            # If e < 0 and term == 0, Nemo will throw a DivideError. 
            # However, for an admissible 6j symbol, denominators never exceed k+1, 
            # so we never divide by Φ_{2N}(z^2) == 0.
            res *= e > 0 ? term^e : inv(term)^abs(e)
        end
    end
    return res
end

# ============================================================
# Exact Cyclo Field Evaluations 
# ============================================================

# Evaluate the triangle coeffs 
function qtricoeff2_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    Δ2 = qtricoeff2_symb(j1, j2, j3, j4, j5, j6)
    return evaluate_cyclofield(Δ2, model)
end

# Evaluate the racah sum arithmetics 
function q6jseries_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    α1 = Int(j1 + j2 + j3); α2 = Int(j1 + j5 + j6) 
    α3 = Int(j2 + j4 + j6); α4 = Int(j3 + j4 + j5)
    β1 = Int(j1 + j2 + j4 + j5); β2 = Int(j1 + j3 + j4 + j6); β3 = Int(j2 + j3 + j5 + j6)
    
    # Removed artificial `model.k` cap. The field arithmetic handles [k+2]_q = 0 exactly.
    zrange = max(α1, α2, α3, α4):min(β1, β2, β3) 
    
    sum_cf = model.K(0)

    @inbounds for z in zrange
        val = q6jsummand_symb(z, α1, α2, α3, α4, β1, β2, β3) 
        sum_cf += evaluate_cyclofield(val, model)
    end
    return sum_cf
end

# Function for quantum6j symbols 
function qracah6j_exact(model::ExactSU2kModel, j1, j2, j3, j4, j5, j6)
    if !qδtet(j1, j2, j3, j4, j5, j6, model.k) 
        return ExactResult(model.k, model.K(0), model.K(0))
    end
    
    Tc2 = qtricoeff2_exact(model, j1, j2, j3, j4, j5, j6)
    Sum_cf = q6jseries_exact(model, j1, j2, j3, j4, j5, j6)
    
    return ExactResult(model.k, Tc2, Sum_cf)
end

# TODO Fulfilled: Standalone Wrapper
"""
    qracah6j_exact(j1, j2, j3, j4, j5, j6, k::Int) -> ExactResult

Creates a temporary ExactSU2kModel and evaluates the exact algebraic 6j symbol.
"""
function qracah6j_exact(j1::Number, j2, j3, j4, j5, j6, k::Int)
    model = ExactSU2kModel(k)
    return qracah6j_exact(model, j1, j2, j3, j4, j5, j6)
end

# ============================================================
# Float Projection (Horner's Method)
# ============================================================

"""
    evaluate_exact(res::ExactResult, prec::Int=256) -> BigFloat

Converts an exact algebraic result into a high-precision float.
"""
function evaluate_exact(res::ExactResult, prec::Int=256)
    # If admissibility failed, it's exactly 0
    if res.pref_sq == 0 && res.sum_cf == 0
        return zero(BigFloat)
    end
    
    setprecision(BigFloat, prec) do
        N = res.k + 2
        # Target z = exp(iπ/N)
        target_z = cispi(big"1.0" / N)
        
        # 1. Evaluate S (Racah Sum) via Horner
        sum_bf = horner_eval(res.sum_cf, target_z)
        
        # 2. Evaluate Δ² (Prefactor squared) via Horner
        pref2_bf = horner_eval(res.pref_sq, target_z)
        
        # 3. Final Assembly
        return real(sqrt(pref2_bf) * sum_bf)
    end
end

"""
    horner_eval(ev::nf_elem, root_val::Complex{BigFloat}) -> Complex{BigFloat}

Evaluates a Nemo number field element at a specific complex root using Horner's scheme.
"""
function horner_eval(ev::nf_elem, root_val::Complex{BigFloat})
    d = degree(parent(ev)) - 1
    d < 0 && return zero(Complex{BigFloat})
    
    res = Complex{BigFloat}(coeff(ev, d))
    for i in (d-1):-1:0
        c = coeff(ev, i)
        # Use Nemo.numerator and Nemo.denominator strictly to avoid Base conflicts
        val = BigFloat(Nemo.numerator(c)) / BigFloat(Nemo.denominator(c))
        res = res * root_val + val
    end
    return res
end