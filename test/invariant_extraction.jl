module InvariantExtractionTests

# Keep interactive Main bindings (for example from a partially included source
# file) from shadowing the package functions exercised by these tests.
using Test, LinearAlgebra, TSSOS, DynamicPolynomials, PermutationGroups, JuMP, COSMO

@testset "Weighted invariant atomic extraction" begin
    fixtures = [
        (:A, [2], 4, [[1],[4]], [1//3,2//3], [(0,),(1,)], :signs, [2,2]),
        (:B, [2,2], 4, [[1,1],[4,1],[1,9]], [1//2,1//3,1//6], [(0,0),(1,0),(0,1)], :signs, [4,4,4]),
        (:C, [1,2,3], 4, [[3,2,0],[4,3,0]], [2//5,3//5], [(0,0,0),(1,0,0)], :permutations, [6,6]),
        (:D, [1,2,3], 5, [[2,0,0],[2,1,0],[3,3,1]], [1//2,1//3,1//6], [(0,0,0),(1,0,0),(0,1,0)], :permutations, [3,3,1]),
        (:E, [2,4,6], 8, [[6,9,4],[14,49,36]], [1//5,4//5], [(0,0,0),(1,0,0)], :signed_permutations, [24,48]),
    ]
    for (id,degrees,t,atoms,weights,B,symmetry,sizes) in fixtures
        @testset "$id" begin
            y = invariant_moments(atoms,weights,degrees,t)
            r = extract_invariant_atoms(y,degrees,t; arithmetic=:exact,basis=B)
            @test r.status == :success
            @test r.rank == length(atoms)
            @test r.diagnostics[:full_representation]
            p = [findfirst(v -> isapprox(v,a;atol=1e-7),atoms) for a in r.atoms]
            @test all(isapprox.(r.atoms,atoms[p];atol=1e-7))
            @test r.weights ≈ weights[p] atol=1e-7
            @test r.basis == B
            @test r.diagnostics[:exact_atoms] == atoms[p]
            @test r.diagnostics[:exact_weights] == weights[p]
            @test r.diagnostics[:exact_full_moments_consistent]
            @test all(A*B == B*A for A in r.multiplication, B in r.multiplication)
            @test all(r.H*T == T'*r.H for T in r.multiplication)
            # Exact multiplication matrices reproduce every supplied rational moment.
            e = [1; zeros(Int,r.rank-1)]
            @test all((r.H * prod(r.multiplication[j]^a[j] for j in eachindex(a)) * e)[1] == value for (a,value) in y)
            lifted = lift_invariant_atoms(r;symmetry,enumerate=true)
            @test getproperty.(lifted,:size) == sizes[p]
            @test all(o -> o.status == :lifted && length(o.points) == o.size,lifted)
            @test all(o -> o.point_weight*o.size ≈ o.weight,lifted)
            @test all(o -> o.label == :candidate && o.feasible === nothing,lifted)
            if id == :A
                @test r.H == [1 3;3 11]
                @test r.multiplication[1] == [0 -4;1 5]
            elseif id == :B
                @test r.multiplication[1] == [0 -4 -1;1 5 1;0 0 1]
                @test r.multiplication[2] == [0 -1 -9;0 1 0;1 1 10]
            elseif id == :C
                @test r.multiplication[1] == [0 -12;1 7]
                @test r.multiplication[2] == r.multiplication[1]-I
                @test iszero(r.multiplication[3])
            elseif id == :D
                @test det(r.H) == 1//36
                @test r.multiplication[3] == r.multiplication[1]-2I
            elseif id == :E
                @test r.multiplication[1] == [0 -84;1 20]
                @test r.multiplication[2] == 5r.multiplication[1]-21I
                @test r.multiplication[3] == 4r.multiplication[1]-20I
            end
            numeric = extract_invariant_atoms(y,degrees,t)
            @test numeric.status == :success
            pnum = [findfirst(v -> isapprox(v,a;atol=1e-6),r.atoms) for a in numeric.atoms]
            @test numeric.weights ≈ r.weights[pnum] atol=1e-6
        end
    end

    @testset "Exponential sign family stays in invariant coordinates" begin
        for n in (3,5,8,10)
            atoms = [ones(Int,n)]
            append!(atoms, [ones(Int,n)+3 .* [j == i for j in 1:n] for i in 1:n])
            weights = [1//2; fill(1//(2n),n)]
            y = invariant_moments(atoms,weights,fill(2,n),4)
            r = extract_invariant_atoms(y,fill(2,n),4)
            @test r.status == :success
            @test r.rank == n+1
            orbits = lift_invariant_atoms(r)
            @test all(o -> o.size == 2^n && o.points === nothing,orbits)
            @test sum(o.size for o in orbits) == (n+1)*2^n
        end
    end

    @testset "Diagnostics and truncation" begin
        y = Dict((0,)=>1//1,(1,)=>1//1,(2,)=>1//1,(3,)=>2//1)
        r = extract_invariant_atoms(y,[2],3)
        @test r.status == :inconsistent_full_moments
        @test r.atoms ≈ [[1.0]]
        @test r.diagnostics[:required_moments_consistent]
        @test !r.diagnostics[:full_moments_consistent]
        @test r.diagnostics[:additional_residual] ≈ 1/3
        @test r.diagnostics[:additional_indices] == [(3,)]
        delete!(y,(3,))
        r = extract_invariant_atoms(y,[2],3)
        @test r.status == :success
        @test !r.diagnostics[:full_degree_coverage]
        @test !r.diagnostics[:full_representation]
        @test_throws ArgumentError extract_invariant_atoms(Dict((0,)=>1),[2],4)
        @test_throws ArgumentError weighted_indices([0],2)
        @test_throws ArgumentError extract_invariant_atoms(y,[2],1)
        @test_throws ArgumentError extract_invariant_atoms(Dict((0,)=>2,(1,)=>2,(2,)=>2),[2],2)
        @test_throws ArgumentError extract_invariant_atoms(Dict((0,)=>1,(1,)=>Inf,(2,)=>1),[2],2)
        bad = Dict((0,)=>1.,(1,)=>0.,(2,)=>-1.)
        @test extract_invariant_atoms(bad,[2],2).status == :not_psd
        nonflat = invariant_moments([[1],[2]],[1//2,1//2],[2],2)
        @test extract_invariant_atoms(nonflat,[2],2).status == :not_flat
        collision = invariant_moments([[1],[4]],[1//3,2//3],[2],4)
        @test extract_invariant_atoms(collision,[2],4;gap_tol=100).status == :spectral_collision
        not_real = extract_invariant_atoms(invariant_moments([[-1]],[1],[2],2),[2],2)
        @test only(lift_invariant_atoms(not_real)).status == :not_realizable
        # No constraint or lower-bound assumptions are manufactured by extraction.
        good = extract_invariant_atoms(collision,[2],4)
        @test all(o -> o.label == :certified_global_minimizers,
            lift_invariant_atoms(good;inequalities=[],equalities=[],objective=x -> (x[1]^2-1)^2*(x[1]^2-4)^2,
                lower_bound=0,lower_bound_certified=true))
        @test all(o -> o.feasible === false,
            lift_invariant_atoms(good;inequalities=[x -> x[1]],equalities=[]))
        custom = lift_invariant_atoms(good;symmetry=:custom,lift=v -> sqrt.(v),invariants=x -> x.^2,
            group=[-1,1],action=(g,x) -> g*x,enumerate=true)
        @test all(o -> o.size == 2 && o.stabilizer == 1 && length(o.points) == 2,custom)
        @test_throws ArgumentError extract_invariant_atoms(collision,[2],4;basis=[(1,),(0,)])
        @test_throws ArgumentError extract_invariant_atoms(collision,[2],4;arithmetic=:unknown)
        nonprob = extract_invariant_atoms(Dict(a=>2v for (a,v) in collision),[2],4;probability=false)
        @test sum(nonprob.weights) ≈ 2
        irrational = extract_invariant_atoms(Dict((0,)=>1//1,(1,)=>0//1,(2,)=>2//1,(3,)=>0//1,(4,)=>4//1),[1],2;arithmetic=:exact)
        @test irrational.status == :success
        @test irrational.diagnostics[:exact_atoms] === nothing
        # Higher precision resolves a separation lost in Float64 moment ranks.
        setprecision(256) do
            for epsilon in (big"1e-3",big"1e-6",big"1e-9")
                atoms = [[big"1",big"1"],[1+epsilon,big"1"]]
                moments = invariant_moments(atoms,[big"0.5",big"0.5"],[2,2],4)
                r = extract_invariant_atoms(moments,[2,2],4)
                @test r.status == :success
                @test r.rank == 2
                @test all(isapprox.(r.atoms,atoms;atol=big"1e-30"))
                double = extract_invariant_atoms(Dict(a=>Float64(v) for (a,v) in moments),[2,2],4)
                @test double.status in (:success,:not_flat,:ill_conditioned_basis,:spectral_collision,:invalid_representation,:not_psd)
                @test haskey(double.diagnostics,:rank_gaps)
            end
        end
    end
end

@testset "Permutation and sign symmetry SDP integration" begin
    @polyvar x[1:2]
    G = PermGroup([perm"(1,2)"])
    f = sum((x.^2 .- 1).^2)
    optimizer = optimizer_with_attributes(COSMO.Optimizer,"eps_abs"=>1e-8,"eps_rel"=>1e-8)
    for TS in (false,"block"), dualize in (false,true)
        opt,data = tssos_symmetry([f,sum(x.^2)-2],x,4,G;numeq=1,signsymmetry=true,TS,dualize,optimizer,QUIET=true)
        @test abs(opt) < 1e-6
        @test length(data.group) == 8
        @test !isempty(data.moments)
        if TS !== false
            # Sparse Gram blocks need not impose PSD on the whole invariant
            # moment block. A dense higher step supplies that condition.
            opt,data = tssos_symmetry(data;TS=false,dualize,optimizer,QUIET=true)
            @test abs(opt) < 1e-6
        end
        r = extract_invariant_atoms(data,[sum(x.^2),prod(x.^2)],x,4;degrees=[2,4],rtol=1e-5,validation_tol=1e-4)
        @test r.status == :success
        @test only(r.atoms) ≈ [2,1] atol=1e-4
        @test only(r.weights) ≈ 1 atol=1e-4
        @test_throws ArgumentError invariant_moments(data,[x[1]],x,4)
    end
    opt,legacy = tssos_symmetry([f],x,2,G;TS=false,optimizer,QUIET=true)
    @test abs(opt) < 1e-6
    opt,data = tssos_symmetry([f;1 .- x.^2],x,2,G;SymmetricConstraint=false,
        signsymmetry=true,semisimple=true,TS="block",optimizer,QUIET=true)
    @test abs(opt) < 1e-6
    @test length(legacy.group) == 2
    # A common sign flip preserves x1*x2; independent flips do not.
    coupled = (sum(x.^2)-2)^2 + (prod(x)-1)^2
    opt,subgroup = tssos_symmetry([coupled],x,2,G;signsymmetry=true,TS=false,optimizer,QUIET=true)
    @test abs(opt) < 1e-6
    @test length(subgroup.group) == 4
    @test TSSOS.normalform(UInt16[1,2],subgroup.group,subgroup.action) !== nothing
    # Equality multipliers, unconsolidated localizers, explicit binary sign columns.
    opt,data = tssos_symmetry([f,1-sum(x.^2)/2,prod(x.^2)-1],x,2,G;
        numeq=1,signsymmetry=Matrix{Int}(I,2,2),TS=false,optimizer,QUIET=true)
    @test abs(opt) < 1e-6
    @test_throws ArgumentError tssos_symmetry([sum(x)],x,1,G;signsymmetry=ones(Int,2,1),optimizer,QUIET=true)
    model = Model(optimizer)
    @variable(model, lower)
    info = add_psatz_symmetry!(model,f-lower,x,[],[],2,G;signsymmetry=true,TS="block",QUIET=true)
    @objective(model,Max,lower)
    optimize!(model)
    @test abs(objective_value(model)) < 1e-6
    @test length(info.group) == 8

    @complex_polyvar z[1:2]
    moduli = z .* conj.(z)
    cf = sum(moduli.^2 + moduli)
    for conjugate in (false,true)
        opt,data = complex_tssos_symmetry([cf],z,2,G;signsymmetry=true,ConjugateBasis=conjugate,TS=false,optimizer,QUIET=true)
        @test abs(opt) < 1e-6
        @test length(data.group) == 8
    end
end

end # module InvariantExtractionTests
