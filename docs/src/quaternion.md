# Quaternion Polynomial Optimization

The `TSSOS.QuaternionSOS` submodule implements the quaternion moment-SOS relaxation,
including correlative and term sparsity. Load its API explicitly:

```julia
using TSSOS.QuaternionSOS
using DynamicPolynomials
using COSMO

@ncpolyvar z[1:2]
pop = [-z[1]*z[2], 1-z[1]*z[2]]
bound = qtssos(pop, z, 1, 1; optimizer=COSMO.Optimizer, QUIET=true)
# bound ≈ -1
```

The variables are ordered as `z[1:n]` followed by their conjugates `z[n+1:2n]`.
Use Hermitian (real-valued) objectives and scalar constraints. `pop[1]` is the
objective, followed by nonnegative constraints and then `numeq`
equalities. The return value is the relaxation bound; `solve=false` computes the
supports without invoking an optimizer and returns `nothing`.

Common keyword arguments:

- `CS="MF"`: correlative sparsity using minimum fill. Also accepts `"MD"`, `"NC"`,
  or `false`. Supply `cliques` to specify variable cliques explicitly.
- `TS=false`: term sparsity; also accepts `"block"`, `"MF"`, or `"MD"`.
- `nb=0`: constrain the first `nb` quaternion variables to unit norm.
- `conjubasis=true`: include conjugate variables in the objective SOS basis.
- `ipart=true`: use quaternion Gram matrices. Set `false` for the real coefficient
  relaxation used in the supplied examples.
- `addrcons=true`: add real commutation identities; `addicons=true` additionally
  adds imaginary commutation identities to their multipliers.
- `rncnumeq=0`, `qncnumeq=0`: counts of trailing noncommutative real and quaternion
  equalities, respectively, within the final `numeq` constraints. Quaternion
  equalities require `ipart=true`. Their multipliers have four independent real
  components, including singleton blocks, and enter the certificate together
  with their quaternion adjoints.
- `optimizer`: a JuMP optimizer factory, or `model`: an empty `JuMP.Model`.
  Without either, TSSOS selects its default optimizer. `QUIET=true` suppresses
  normal progress output.

`fsupp` and `fcoe` may supply an explicit quaternion objective when `ipart=true`.
Include the complete Hermitian representation, with adjoint terms and their
coefficients. Each support term is a vector of three `Vector{UInt16}` parts: paired variable
indices for scalar norm factors, their conjugate indices, and the remaining
ordered word. Supply both arrays with matching lengths. This representation is
useful when quaternion coefficients occur between noncommuting variables.

`quaternion_to_real(pop, z)` returns real polynomials and `4n` commuting variables,
in blocks containing the real, i, j, and k coordinates. It takes the real part
after substitution. This permits comparison with `tssos` or `cs_tssos`; see
`example/quaternion.jl`.

```@docs
TSSOS.QuaternionSOS.qtssos
TSSOS.QuaternionSOS.quaternion_to_real
```
