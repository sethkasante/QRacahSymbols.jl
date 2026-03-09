# Theory & Architecture

The **QRacahSymbols.jl** package implements the $q$-deformed Racah formula for the $SU(2)_k$ quantum group.


The core of the library is the evaluation of the quantum $6j$-symbol 
$$\begin{Bmatrix} j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix}_q,$$ 
which represents the associativity of the fusion of three anyons in a topological quantum field theory.



## The Kirillov-Reshetikhin Formula
We evaluate the quantum Racah $W$-coefficient ($6j$-symbol) using the following alternating sum:

$$\;\;\; \quad \begin{Bmatrix} j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix}_q = \Delta(j_1, j_2, j_3) \Delta(j_1, j_5, j_6) \Delta(j_2, j_4, j_6) \Delta(j_3, j_4, j_5) \sum_z \frac{(-1)^z [z+1]_q!}{\prod_i [z - \alpha_i]_q! \prod_j [\beta_j - z]_q!}$$

Where:
* $[n]_q = \frac{q^{n/2} - q^{-n/2}}{q^{1/2} - q^{-1/2}}, \;\;\;\;\;\;\; \text{with} \;\;\;\; q = e^{i \frac{2\pi}{k+2}}$ is the quantum integer.
* $\Delta(a, b, c) = \left( \frac{[a+b-c]_q! \; [a-b+c]_q! \; [-a+b+c]_q!}{[a+b+c+1]_q!}\right)^\frac12$ is the quantum triangle coefficient.
* The sum is stabilized in the `Numeric` engine using a **log-sum-exp** shift to prevent floating-point overflow.

## The Catastrophic Cancellation Problem
Topological Quantum Field Theories (TQFTs) rely heavily on $\text{SU(2)}_k$ modular tensor category data. The core invariant is given by the quantum Racah $W$-coefficient ($6j$-symbol), an alternating sum of massive quantum factorials. 

In standard floating-point arithmetic (`Float64`), this alternating series suffers from **catastrophic cancellation** at high spins. For example, at $j=500$, the terms in the sum can reach $10^{60}$, completely obliterating the precision of the final answer (which is often on the order of $10^{-5}$).

## The CycloMonomial Solution
To bypass this limit, `QRacahSymbols.jl` implements a `CycloMonomial` architecture. Because the quantum deformation parameter $q$ is a root of unity, the quantum integers $[n]_q$ can be factored exactly into **cyclotomic polynomials** $\Phi_d(q)$.

Instead of computing massive intermediate floating-point numbers, the Generic Engine represents the Racah sum purely as arrays of prime-power exponents. The massive factorials are cancelled algebraically $O(1)$ *before* any numerical evaluation takes place.

### Cyclotomic Factorization
In the `Generic` mode, we exploit the fact that quantum integers $[n]_q$ are products of *cyclotomic polynomials* $\Phi_d(q)$. Thus, any quantum factorial can be represented as cyclotomic monomials:

$$[n]_q! =  s \cdot z^{p} \cdot \prod_d \Phi_d(q)^{e_d}$$

with integer exponents, where $z = q^{1/2}$. This allows for exact symbolic manipulation without ever choosing a specific value for $q$. The quantum $\{6j\}_q$-symbols is therefore a linear combination of products of cyclotomic monomials. 



## Topological Category Data

Beyond the $\{6j\}_q$-symbols, the package natively evaluates the full suite of TQFT structural constants required for state-sum models: 
* Quantum $\{3j\}_q$-symbols 
* Quantum Dimensions ($[2j+1]_q$)
* F-symbols (Unitary fusion matrices)
* R-matrices (Braiding phases)
* G-symbols (Tetrahedral weights for Turaev-Viro models)
