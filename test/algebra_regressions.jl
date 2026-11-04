using Test
using LinearAlgebra
using COSMO
using DynamicPolynomials
using JuMP
using TSSOS

@testset "Algebraic regression cases" begin
    @testset "Cross-clique unit-modulus constraints" begin
        @complex_polyvar z[1:2]
        squares = z .* conj.(z)
        g = prod(squares) - 1
        for numeq in (0, 1), TS in ("block", false)
            opts = (; cliques=[[1],[2]], normality=0, numeq, TS, QUIET=true)
            bound, _, _ = complex_cs_tssos([sum(squares), g], z, 2; nb=2, opts...)
            @test bound ≈ 2 atol=1e-6
            # On |z₁|=1 the constraint still imposes |z₂|² ≥ 1 (or = 1).
            bound, _, _ = complex_cs_tssos(Poly{Float64}[squares[2], g], z, 2; nb=1, opts...)
            @test bound ≈ 1 atol=1e-6
        end
    end

    @testset "Matrix weights with separated duplicates" begin
        left = [0.3 0.1; 0.1 0.2]
        right = [0.1 0.0; 0.0 0.4]
        mass, first_moment = left + right, right - left
        moment = [mass first_moment; first_moment mass]
        for (atoms, expected) in (([[-1.], [1.], [-1.]], [left, right]),
                                  ([[1.], [-1.], [1.00001]], [right, left]))
            weights = extract_weight_matrix(moment, 1, 1, 2, atoms)
            @test length(weights) == 2
            @test all(isapprox.(weights, expected))
        end
        @test_throws ArgumentError extract_weight_matrix(moment, 1, 1, 2, [[-1.], [0.], [1.]])
    end

    @testset "Dense sum-of-ratios optimizer" begin
        @polyvar x
        q = [1.0+x^2, 1.0+x^2]
        @test SumOfRatios(q, q, [], [], [x], 1; QUIET=true) ≈ 2 atol=1e-6
        solver = optimizer_with_attributes(COSMO.Optimizer, "eps_abs"=>1e-8, "eps_rel"=>1e-8)
        @test SumOfRatios(q, q, [], [], [x], 1; optimizer=solver, dualize=true, QUIET=true) ≈ 2 atol=1e-6
    end

    @testset "Polynomial identities" begin
        @polyvar x
        @complex_polyvar z
        for (builder, variable, target) in (
            (add_psatz!, x, -x^4),
            (add_complex_psatz!, z, -z^2*conj(z)^2),
        )
            model = Model(COSMO.Optimizer)
            set_silent(model)
            builder(model, convert(Poly{Float64}, target), [variable], [], [], 1; TS=false, QUIET=true)
            optimize!(model)
            @test termination_status(model) == MOI.INFEASIBLE
        end

        model = Model(COSMO.Optimizer)
        set_silent(model)
        @variable(model, t)
        add_psatz!(model, 1.0 + t*x^4, [x], [], [], 1; TS=false, QUIET=true)
        @objective(model, Max, t)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(t) ≈ 0 atol=1e-6

        model = Model(COSMO.Optimizer)
        set_silent(model)
        @variable(model, scale)
        add_psatz_cheby!(model, 1.0 - (1e-9*scale)*x^2, [x], [], [], 1; TS=false, QUIET=true)
        # Check the affine coefficient before solving to distinguish a missing
        # coefficient equation from a solver tolerance issue.
        rows = all_constraints(model, AffExpr, MOI.EqualTo{Float64})
        @test any(c -> !iszero(normalized_coefficient(c, scale)), rows)
        fix(scale, 1e9)
        optimize!(model)
        @test termination_status(model) == MOI.INFEASIBLE
    end

    @testset "Quotient algebra" begin
        @polyvar x y
        f = x^2 + (y-2)^2
        expected = 1 + (y-2)^2
        @test TSSOS.binary_normalform(f, [x,y], 1) == expected
        for input in ([f], [poly(f, [x,y])])
            variables = input[1] isa Poly ? [x,y] : 2
            opt, _, _ = cs_tssos(input, variables, 1; nb=1, CS=false, QUIET=true)
            @test opt ≈ 1 atol=1e-6
        end
        for (pop, options) in (([f], (; nb=1)),
                               ([f, x^2-1], (; numeq=1)),
                               ([f, x^2-1], (; numeq=1, GroebnerBasis=true)))
            opt, _, data = tssos(pop, [x,y], 1; options..., QUIET=true)
            @test opt ≈ 1 atol=1e-6
            @test data.pop[1] == f
        end

        opt, _, _ = tssos([y^2, x-y], [x,y], 1; nb=1, numeq=1, GroebnerBasis=true, QUIET=true)
        @test opt ≈ 1 atol=1e-6
        opt, _, _ = tssos([f, x^2-1], [x,y], 1; nb=1, numeq=1, TS="signsymmetry", QUIET=true)
        @test opt ≈ 1 atol=1e-6

        for (relation, expected_square) in ((x^2-1, 1.0), (x^2-x-1, 1.25))
            info = add_psatz!(Model(), convert(Poly{Float64}, 1+x^2), [x], [], [relation], 1;
                              TS=false, GroebnerBasis=true, QUIET=true)
            first_moment = relation == x^2-1 ? 0.0 : 0.25
            moments = [isempty(mon) ? 1.0 : first_moment for mon in info.tsupp]
            @test get_moment_matrix(moments, info)[1] == [1 first_moment; first_moment expected_square]
        end
    end

    @testset "Complex pruning and first moments" begin
        @complex_polyvar z[1:2]
        for shift in (0, 3)
            f = shift + sum(v*conj(v) for v in z)
            for prune in (false, true)
                opt, _, _ = complex_cs_tssos(Poly{Float64}[f], z, 2; cliques=[[1],[2]],
                                              normality=0, reducebasis=prune, QUIET=true)
                @test opt ≈ shift atol=1e-6
            end
        end
        for normality in (0, 1)
            model = Model(COSMO.Optimizer)
            set_silent(model)
            f = convert(Poly{Float64}, z[1]*conj(z[1]))
            add_complex_psatz!(model, f, z[1:1], [], [], 1; ConjugateBasis=true, normality, QUIET=true)
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
        end
        lower, _, _ = complex_tssos(Poly{Float64}[z[1]+conj(z[1])], z[1:1], 1;
                                      nb=1, normality=0, reducebasis=true, QUIET=true)
        @test lower ≈ -2 atol=1e-6

        a = 1 + 2im
        one_mon, z_mon = UInt16[], UInt16[1]
        basis = [one_mon, z_mon]
        supports = [(one_mon, one_mon), (one_mon, z_mon), (z_mon, z_mon)]
        # The off-diagonal dual stores twice the real and imaginary first moments.
        exact = [TSSOS.complex_moment((u,v), [1,2,5], [4], supports, [supports[2]]) for u in basis, v in basis]
        @test exact == [1 conj(a); a abs2(a)]

        f = (z[1]-a)*(conj(z[1])-conj(a))
        opt, sol, data = complex_cs_tssos([cpoly(f, z[1:1])], 1, 1; CS=false, normality=0,
                                          solution=true, QUIET=true)
        @test opt ≈ 0 atol=1e-6
        @test data.flag == 0
        # The objective error is quadratic in the distance to the minimizer.
        @test only(only(sol)) ≈ a atol=1e-4
        @test data.moment[1][1][2,1] ≈ a atol=1e-4
    end

    @testset "Matrix coefficients and moments" begin
        @polyvar x
        P = typeof(1.0+x)
        F = [P[3 4; 4 3], P[1 0; 0 1]]
        G = [P[1 1; 1 1]]
        @test sparseobj([1.0], F, G, [x], 0; TS=false, QUIET=true)[1] ≈ 1 atol=1e-6

        A = [1 2; 2 5; 3 1; 1 4]
        Q = A*A'
        constraint = P[x^2 1+x; 1+x 2]
        for j = 1:2, k = j:2
            block = Q[2j-1:2j, 2k-1:2k]
            @test TSSOS.multiplier_contraction(triu(Q), constraint, j, k) == sum(block .* constraint)
        end

        scalar_pmi(g) = [reshape(P[-1], 1, 1), reshape(P[g], 1, 1)]
        opt, data = LinearPMI([1.0], scalar_pmi(1-x^4), P[], [x], 1; TS=false, QUIET=true)
        @test data.SDP_status == MOI.INFEASIBLE
        @test isnothing(opt)
        @test isnothing(data.moment)
        opt, data = LinearPMI([1.0], scalar_pmi(1+x^2), P[], [x], 1; TS=false, QUIET=true)
        @test opt ≈ 1 atol=1e-6
        @test only(data.moment) ≈ [1 0; 0 0] atol=1e-6

        opt, _, data = tssos(P[1+x^2 0; 0 2+2x^2], [1.0-x^2], [x], 1; Moment=true, QUIET=true)
        @test opt ≈ 1 atol=1e-6
        @test size(only(data.moment)) == (4,4)
        @test only(data.moment)[1,2] == 0
        basis = [[UInt16[], UInt16[1]]]
        supports = [[UInt16[], UInt16[1,1]], [UInt16[1]], [UInt16[], UInt16[1], UInt16[1,1]]]
        M = only(get_mmoment([[2.,4.], [6.], [7.,8.,9.]], supports, 1, basis, 2))
        @test M == [2 0 0 3; 0 7 3 8; 0 3 4 0; 3 8 0 9]
    end

    @testset "Certification on the original domain" begin
        @polyvar x
        f = Poly{Float64}[x^2]
        @test !TSSOS.check_optimality([0.0], 0.0, f, [x]; nb=1)
        @test !TSSOS.check_optimality([NaN], 1.0, f, [x]; nb=1)
        @test isnothing(extract_solutions(Diagonal([1., 0.]), 1; x=[x], pop=f, lb=0., check=true, nb=1))
        @test !TSSOS.check_optimality([1e200], 0., Poly{Float64}[0*x, x^2], [x])
        opt, sol, data = tssos(Poly{Float64}[-x^2], [x], 1; solution=true, QUIET=true)
        @test isnothing(opt) && isnothing(sol) && isnothing(data.moment)
        @test data.SDP_status == MOI.INFEASIBLE
        opt, sol, data = cs_tssos(f, [x], 1; CS=false, nb=1, solution=true, QUIET=true)
        @test opt ≈ 1 atol=1e-6
        @test data.flag == 0
        @test abs(only(sol)) ≈ 1 atol=1e-6

        @complex_polyvar z[1:2]
        opt, sol, data = complex_tssos(Poly{Float64}[z[1]*conj(z[1])], z[1:1], 1;
                                       nb=1, normality=0, solution=true, QUIET=true)
        @test opt ≈ 1 atol=1e-6
        @test data.flag == 0
        @test abs(only(sol)) ≈ 1 atol=1e-6
        # The first sparse relaxation has a known gap. A locally solved
        # refinement may still return a useful feasible candidate.
        a,b = z
        pop = [3-a*conj(a)-0.5im*a*conj(b)^2+0.5im*b^2*conj(a), b+conj(b),
               -1+a*conj(a)-0.25*a^2-0.25*conj(a)^2, -3+a*conj(a)+b*conj(b), im*(b-conj(b))]
        lower, _, data = complex_tssos(pop, z, 2; numeq=3, solution=true, QUIET=true)
        @test lower ≈ -0.4165019 atol=1e-5
        @test data.flag == 1
        lower, candidate, data = complex_tssos(data; solution=true, QUIET=true)
        @test lower ≈ 0.4281746 atol=1e-5
        certified = !isnothing(candidate) && TSSOS.check_optimality(candidate, lower, pop, z; numeq=3)
        @test (data.flag == 0) == certified
    end

    @testset "Binary multiplication" begin
        @polyvar x[1:2]
        for n in (1, 2)
            d = n+1
            atoms = collect(Iterators.product(fill((-1,1), n)...))
            basis = get_basis(n, d)
            V = [prod((atom[i] for i in mon); init=1) for mon in basis, atom in vec(atoms)]
            moment = Symmetric(V*V' / length(atoms))
            solutions = extract_solutions(moment, d; x=x[1:n], nb=n)
            @test length(solutions) == length(atoms)
            @test all(atom -> any(s -> isapprox(s, collect(atom); atol=1e-6), solutions), atoms)
        end
    end

    @testset "Christoffel polynomials" begin
        @polyvar x[1:2]
        M = Diagonal([1.0,1e-4])
        _, polynomial, rank, _, _ = construct_CDK(x[1:1], 1, M; threshold=1e-6)
        @test rank == 2
        @test polynomial == 1+10000*x[1]^2
        @test construct_CDK(x[1:1], 1, M)[3] == 1
        kernel, polys, ranks, minima, nullities = construct_CDK_cs(x, 1, [M, Diagonal([1., 0.5])], [[1],[2]])
        @test only(kernel[1]) == x[1]^2
        @test isempty(kernel[2])
        @test polys == [1, 1+2*x[2]^2]
        @test ranks == [1,2]
        @test minima == [1e-4,0.5]
        @test nullities == [1,0]
        @test ranks !== nullities
    end
end
