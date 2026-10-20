using PermutationGroups
using LinearAlgebra

# Use the solver from issue #15 for dense PSD blocks. Clarabel 0.11.1's
# symmetric Kronecker products still contain a Float64 sqrt(2) constant.
function high_precision_sdp_model()
    model = GenericModel{BigFloat}(ClusteredLowRankSolver.Optimizer)
    set_optimizer_attribute(model, "duality_gap_threshold", big"1e-30")
    return model
end

@testset "High-precision solver interfaces" begin
    setprecision(BigFloat, 256) do
        @polyvar x[1:2]
        c, a = 1 + big"1e-20", big"0.25"
        @testset "tssos" begin
            for TS in ("block", false)
                opt, sol, data = tssos([x[1]^2+c, 1-x[1]^2, x[1]-a], x[1:1], 1;
                    numeq=1, TS, model=high_precision_model(), Gram=true, solution=true, QUIET=true)
                @test opt isa BigFloat
                @test abs(opt-c-a^2) < big"1e-25"
                @test all(eltype(m) == BigFloat for m in data.moment)
                @test all(eltype(g) == BigFloat for blocks in data.GramMat for g in blocks)
                @test all(eltype(m) == BigFloat for m in data.multiplier)
                @test sol !== nothing
            end
            # The Newton-polytope LP and final SDP share the supplied optimizer.
            opt, _, data = tssos(x[1]^4+c, x[1:1]; model=high_precision_model(), Gram=true, QUIET=true)
            @test abs(opt-c) < big"1e-25"
            @test eltype(data.moment[1]) == BigFloat
            opt, _, data = tssos([x[1]^2+c, 1-x[1]^2], x[1:1], 1;
                model=high_precision_model(), normality=true, Gram=true, QUIET=true)
            @test abs(opt-c) < big"1e-25"
            opt, sol, data = tssos(data; TS=false, model=high_precision_model(), Gram=true, solution=true, QUIET=true)
            @test opt isa BigFloat
            @test abs(opt-c) < big"1e-25"
            @test sol !== nothing
            for options in ((;nb=1), (;numeq=1, GroebnerBasis=true))
                pop = haskey(options, :nb) ? [x[1]^2+c] : [x[1]^2+c, x[1]-a]
                opt, _, data = tssos(pop, x[1:1], 1; options..., model=high_precision_model(), QUIET=true)
                @test abs(opt-c-(haskey(options, :nb) ? 1 : a^2)) < big"1e-25"
            end
            opt, _, _ = tssos(x[1]^2+c, x[1:1]; feasibility=true,
                model=high_precision_model(), QUIET=true)
            @test abs(opt-(1+c)) < big"1e-25"
        end

        @complex_polyvar z[1:2]
        b = a + big"0.125"*im
        @testset "complex_tssos / complex_cs_tssos" begin
            for solver in (complex_tssos, complex_cs_tssos), TS in ("block", false), ConjugateBasis in (false, true)
                n = solver === complex_tssos ? 1 : 2
                vars = z[1:n]
                pop = [sum((zi-b)*conj(zi-b) for zi in vars)+c; 1 .- vars.*conj.(vars)]
                opt, sol, data = solver(pop, vars, 1; TS, ConjugateBasis,
                    model=high_precision_sdp_model(), Gram=true, solution=!ConjugateBasis, QUIET=true)
                @test opt isa BigFloat
                @test abs(opt-c) < big"1e-25"
                @test all(eltype(m) == Complex{BigFloat} for clique in data.moment for m in clique)
                @test all(eltype(g) == Complex{BigFloat} for clique in data.GramMat for blocks in clique for g in blocks)
                if !ConjugateBasis
                    @test sol !== nothing
                end
            end
            pop = [z[1]*conj(z[1])+c, 1-z[1]*conj(z[1])]
            opt, _, data = complex_tssos(pop, z[1:1], 1; normality=0,
                model=high_precision_sdp_model(), Gram=true, QUIET=true)
            @test abs(opt-c) < big"1e-25"
            @test eltype(data.moment[1][1]) == BigFloat
            opt, _, data = complex_tssos(data; TS=false, model=high_precision_sdp_model(), Gram=true, QUIET=true)
            @test opt isa BigFloat
            @test abs(opt-c) < big"1e-25"
            pop = [(z[1]-b)*conj(z[1]-b)+c, 1-z[1]*conj(z[1]), z[1]+conj(z[1])-2a]
            opt, _, data = complex_cs_tssos(pop, z[1:1], 1; numeq=1,
                model=high_precision_sdp_model(), Gram=true, QUIET=true)
            @test abs(opt-c) < big"1e-25"
            @test eltype(coefficients(data.multiplier[1][1])) == Complex{BigFloat}
        end

        G = PermGroup([perm"(1,2)"])
        @testset "symmetry and certificates" begin
            pop = [sum((xi-a)^2 for xi in x)+c, 1-sum(x.^2)]
            for TS in ("block", false), semisimple in (false, true)
                opt, data = tssos_symmetry(pop, x, 1, G; TS, semisimple,
                    model=high_precision_sdp_model(), QUIET=true)
                @test opt isa BigFloat
                @test abs(opt-c) < big"1e-25"
                @test all(eltype(p.coe) == BigFloat for bas in data.basis[1] for p in bas)
                @test all(abs(sum(abs2, p.coe)-1) < big"1e-60" for bas in data.basis[1] for p in bas)
                @test all(g isa BigFloat || eltype(g) == BigFloat for gs in data.GramMat for bs in gs for g in bs)
            end
            for ConjugateBasis in (false, true), shift in (a, b)
                pop = [sum((zi-shift)*conj(zi-shift) for zi in z)+c, 1-sum(z.*conj.(z))]
                opt, data = complex_tssos_symmetry(pop, z, 1, G; ConjugateBasis,
                    model=high_precision_sdp_model(), QUIET=true)
                @test opt isa BigFloat
                @test abs(opt-c) < big"1e-25"
                @test all(eltype(p.coe) in (BigFloat, Complex{BigFloat}) for bas in data.basis[1] for p in bas)
            end
            model = high_precision_sdp_model()
            set_silent(model)
            @variable(model, lower)
            info = add_psatz_symmetry!(model, sum((xi-a)^2 for xi in x)+c-lower,
                x, [1-sum(x.^2)], [], 1, G; QUIET=true)
            @objective(model, Max, lower)
            optimize!(model)
            @test abs(objective_value(model)-c) < big"1e-25"
            @test info.coe_type == BigFloat
            @test all(v isa GenericVariableRef{BigFloat} for v in all_variables(model))

            # Non-abelian minimal projections and non-real characters exercise
            # both branches of the high-precision symmetry decomposition.
            @polyvar y[1:3]
            @complex_polyvar w[1:3]
            for group in (PermGroup([perm"(1,2)", perm"(1,2,3)"]), PermGroup([perm"(1,2,3)"]))
                opt, data = tssos_symmetry([sum((yi-a)^2 for yi in y)+c, 1-sum(y.^2)], y, 1, group;
                    model=high_precision_sdp_model(), QUIET=true)
                @test abs(opt-c) < big"1e-25"
                @test all(abs(sum(abs2, p.coe)-1) < big"1e-60" for bas in data.basis[1] for p in bas)
                opt, data = complex_tssos_symmetry([sum((wi-b)*conj(wi-b) for wi in w)+c,
                    1-sum(w.*conj.(w))], w, 1, group; model=high_precision_sdp_model(), QUIET=true)
                @test abs(opt-c) < big"1e-25"
                @test all(abs(sum(abs2, p.coe)-1) < big"1e-60" for bas in data.basis[1] for p in bas)
            end
            for solver in (tssos_symmetry, complex_tssos_symmetry)
                vars = solver === tssos_symmetry ? x : z
                f = solver === tssos_symmetry ? sum(vars.^2)+c : sum(vars.*conj.(vars))+c
                opt, data = solver([f], vars, 1, G; model=high_precision_sdp_model(), QUIET=true)
                opt, data = solver(data; TS=false, model=high_precision_sdp_model(), QUIET=true)
                @test opt isa BigFloat
                @test abs(opt-c) < big"1e-25"
            end
            opt, data = tssos_symmetry([sum(x.^2)+c; x.^2 .- a], x, 1, G; numeq=2,
                SymmetricConstraint=false, model=high_precision_sdp_model(), QUIET=true)
            @test abs(opt-c-2a) < big"1e-25"
            @test all(eltype(m) == BigFloat for m in data.multiplier)
        end
    end
end

@testset "Float64 compatibility of extended interfaces" begin
    @polyvar x[1:2]
    @complex_polyvar z[1:2]
    G = PermGroup([perm"(1,2)"])
    for solver in (tssos, complex_tssos, complex_cs_tssos, tssos_symmetry, complex_tssos_symmetry)
        vars = solver in (tssos, tssos_symmetry) ? x : z
        f = vars === x ? sum((xi-0.25)^2 for xi in x) : sum((zi-0.25)*conj(zi-0.25) for zi in z)
        pop = [f, vars === x ? 1-sum(x.^2) : 1-sum(z.*conj.(z))]
        args = solver in (tssos_symmetry, complex_tssos_symmetry) ? (pop, vars, 1, G) : (pop, vars, 1)
        result = solver(args...; optimizer=Clarabel.Optimizer, QUIET=true)
        @test result[1] isa Float64
        @test abs(result[1]) < 1e-6
    end
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, lower)
    add_psatz_symmetry!(model, sum((xi-0.25)^2 for xi in x)-lower, x, [1-sum(x.^2)], [], 1, G; QUIET=true)
    @objective(model, Max, lower)
    optimize!(model)
    @test abs(objective_value(model)) < 1e-6
end
