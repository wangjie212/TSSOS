using DynamicPolynomials
# import DynamicPolynomials as DP

# # Defining action on polynomials by acting on terms and monomials:
# function action(a::Action, el::GroupElement, poly::DP.AbstractPolynomial)
#     return sum(action(a, el, term) for term in DP.terms(poly))
# end

# function action(a::Action, el::GroupElement, term::DP.AbstractTerm)
#     return DP.coefficient(term) * action(a, el, DP.monomial(term))
# end

struct VariablePermutation{V} <: ByPermutations
    variables::V
end

function action(
    a::VariablePermutation,
    g::AP.AbstractPermutation,
    m::AbstractMonomial,
)
    v = a.variables
    return m(v => action(a, g, v))
end

# this is a general linear action that can be induced
# from the action on monomials
# abstract type OnMonomials <: ByLinearTransformation end

# function decompose(
#     k::DP.AbstractPolynomialLike,
#     hom::InducedActionHomomorphism,
# )
#     # correct only if basis(hom) == monomials
#     I = _int_type(hom)
#     indcs = I[hom[mono] for mono in DP.monomials(k)]
#     coeffs = DP.coefficients(k)

#     return indcs, coeffs
# end
