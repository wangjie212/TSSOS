using Test
using JuMP
using DynamicPolynomials
using Clarabel
using ClusteredLowRankSolver
using TSSOS

function high_precision_model()
    model = GenericModel{BigFloat}(Clarabel.Optimizer{BigFloat})
    for name in ("tol_gap_abs", "tol_gap_rel", "tol_feas")
        set_optimizer_attribute(model, name, big"1e-30")
    end
    # Clarabel's regularization/refinement defaults target Float64 accuracy.
    for name in ("static_regularization_constant", "dynamic_regularization_delta",
                 "iterative_refinement_reltol", "iterative_refinement_abstol")
        set_optimizer_attribute(model, name, big"1e-40")
    end
    set_optimizer_attribute(model, "dynamic_regularization_eps", big"1e-50")
    return model
end

@testset "High-precision cs_tssos (issue #15)" begin
    setprecision(BigFloat, 256) do
        @polyvar x[1:1]
        # This offset is lost if any coefficient passes through Float64.
        c = 1 + big"1e-20"
        a = big"0.25"
        pop = [x[1]^2 + c, 1 - x[1]^2, x[1] - a]
        for TS in ("block", false)
            model = high_precision_model()
            opt, sol, data = cs_tssos(pop, x, 1; model, numeq=1,
                TS, Gram=true, solution=true, QUIET=true)
            # Equality constraints make this SDP degenerate; check accuracy even
            # when Clarabel reports reduced-accuracy termination.
            @test data.SDP_status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
            @test opt isa BigFloat
            @test abs(opt - (c + a^2)) < big"1e-25"
            @test all(eltype(m) == BigFloat for m in data.moment)
            @test all(eltype(g) == BigFloat for clique in data.GramMat for blocks in clique for g in blocks)
            @test all(eltype(v) == BigFloat for clique in data.multiplier for v in clique)
            @test sol !== nothing
            if sol !== nothing
                extracted = TS === false ? sol[1] : sol
                @test extracted[1] ≈ Float64(a) atol=1e-6
            end
            @test data.flag == 0
            # Extraction remains a Float64 heuristic, without changing stored moments.
            @test eltype(TSSOS.extraction_moments(data.moment)[1]) == Float64
            @test eltype(data.moment[1]) == BigFloat
            @test all(v isa GenericVariableRef{BigFloat} for v in all_variables(model))
            @test constraint_object(model[:con]).func[1].constant == -c
        end

        # Exercise the solver mentioned in the report as well as Clarabel.
        clrs_model = GenericModel{BigFloat}(ClusteredLowRankSolver.Optimizer)
        set_optimizer_attribute(clrs_model, "duality_gap_threshold", big"1e-30")
        opt, sol, data = cs_tssos([(x[1] - a)^2 + c, 1 - x[1]^2], x, 1;
            model=clrs_model, Gram=true, solution=true, QUIET=true)
        @test opt isa BigFloat
        @test abs(opt - c) < big"1e-25"
        @test eltype(data.moment[1]) == BigFloat
        @test sol !== nothing

        # Integer input, support-form input, and higher TS steps use the supplied model too.
        integer_pop = [x[1]^2, 1 - x[1]^2]
        npop = [poly(p, x) for p in integer_pop]
        opt, _, data = cs_tssos(npop, 1, 1; model=high_precision_model(),
            Gram=true, MomentOne=true, QUIET=true)
        @test opt isa BigFloat
        @test abs(opt) < big"1e-25"
        @test eltype(data.moment[1]) == BigFloat
        old_blocks = deepcopy(data.blocks)
        opt, _, data = cs_tssos(data; model=high_precision_model(),
            TS=false, Gram=true, QUIET=true)
        @test data.blocks != old_blocks
        @test opt isa BigFloat
        @test abs(opt) < big"1e-25"
        @test eltype(data.moment[1]) == BigFloat

        # Moment reconstruction must retain digits beyond Float64, without rounding.
        dual_values = BigFloat[1, a + big"1e-30", a^2]
        supports = [UInt16[], UInt16[1], UInt16[1, 1]]
        moments = get_moment(dual_values, supports, [UInt16[1]], 1, [1])
        @test moments[1][1, 2] == dual_values[2]
        @test moments[1][1, 2] != BigFloat(Float64(dual_values[2]))

        # No solver is needed for a structure-only call, even with solution=true.
        opt, sol, _ = cs_tssos(integer_pop, x, 1; solve=false, solution=true, QUIET=true)
        @test opt === nothing
        @test sol === nothing

        @test_throws ArgumentError cs_tssos(integer_pop, x, 1;
            model=high_precision_model(), optimizer=Clarabel.Optimizer, QUIET=true)
        @test_throws ArgumentError cs_tssos(integer_pop, x, 1;
            model=high_precision_model(), dualize=true, QUIET=true)
        nonempty = high_precision_model()
        @variable(nonempty, existing)
        @test_throws ArgumentError cs_tssos(integer_pop, x, 1; model=nonempty, QUIET=true)

        # Existing optimizer factories and supplied Float64 models remain supported.
        for options in ((; optimizer=Clarabel.Optimizer),
                        (; optimizer=Clarabel.Optimizer, dualize=true),
                        (; model=Model(Clarabel.Optimizer)))
            opt, _, data = cs_tssos(integer_pop, x, 1; options..., Gram=true, QUIET=true)
            @test opt isa Float64
            @test abs(opt) < 1e-6
            @test eltype(data.moment[1]) == Float64
        end
    end
end
