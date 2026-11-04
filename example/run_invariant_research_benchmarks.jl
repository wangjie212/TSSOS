module InvariantResearchBenchmarks

using TSSOS, LinearAlgebra, Statistics, Random, TOML, Dates, SHA

BLAS.set_num_threads(1)
const SEED = 0
const OUT = abspath(isempty(ARGS) ? joinpath(@__DIR__,"..","output","symmetry_benchmarks") : ARGS[1])
mkpath(OUT)

function fixture(id; n=3)
    id == "A" && return (;id,n=1,degrees=[2],t=4,atoms=[[1],[4]],weights=[1//3,2//3],B=[(0,),(1,)],symmetry=:signs)
    id == "B" && return (;id,n=2,degrees=[2,2],t=4,atoms=[[1,1],[4,1],[1,9]],weights=[1//2,1//3,1//6],B=[(0,0),(1,0),(0,1)],symmetry=:signs)
    id == "C" && return (;id,n=3,degrees=[1,2,3],t=4,atoms=[[3,2,0],[4,3,0]],weights=[2//5,3//5],B=[(0,0,0),(1,0,0)],symmetry=:permutations)
    id == "D" && return (;id,n=3,degrees=[1,2,3],t=5,atoms=[[2,0,0],[2,1,0],[3,3,1]],weights=[1//2,1//3,1//6],B=[(0,0,0),(1,0,0),(0,1,0)],symmetry=:permutations)
    id == "E" && return (;id,n=3,degrees=[2,4,6],t=8,atoms=[[6,9,4],[14,49,36]],weights=[1//5,4//5],B=[(0,0,0),(1,0,0)],symmetry=:signed_permutations)
    id == "F" || error("unknown fixture")
    atoms = [ones(Int,n)]
    append!(atoms,[ones(Int,n)+3 .* [i==j for j in 1:n] for i in 1:n])
    B = [ntuple(_ -> 0,n); [ntuple(j -> Int(i==j),n) for i in 1:n]]
    return (;id,n,degrees=fill(2,n),t=4,atoms,weights=[1//2;fill(1//(2n),n)],B,symmetry=:signs)
end

function objective(f,x)
    a = x.^2
    f.id == "A" && return ((a[1]-1)*(a[1]-4))^2
    f.id == "B" && return ((a[1]-1)*(a[1]-4))^2+((a[2]-1)*(a[2]-9))^2+((a[1]-1)*(a[2]-1))^2
    if f.id in ("C","D","E")
        e = TSSOS._elementary(f.id == "E" ? a : x)
        f.id == "C" && return (e[2]-e[1]+1)^2+e[3]^2+((e[1]-3)*(e[1]-4))^2
        f.id == "D" && return ((e[1]-2)*(e[1]-3))^2+((e[1]-2)*(e[2]-3))^2+((e[1]-3)*e[2]*(e[2]-1))^2+(e[3]-e[1]+2)^2
        return ((e[1]-6)*(e[1]-14))^2+(e[2]-5e[1]+21)^2+(e[3]-4e[1]+20)^2
    end
    return sum(((v-1)*(v-4))^2 for v in a)+sum(((a[i]-1)*(a[j]-1))^2 for i in 1:f.n for j in i+1:f.n)
end

function measure(f; repeats=31)
    f(); f()
    runs = [@timed(f()) for _ in 1:repeats]
    times = getproperty.(runs,:time)
    value = runs[end].value
    stages = hasproperty(value,:diagnostics) && haskey(value.diagnostics,:timings_seconds) ?
        Dict(string(k)=>median([run.value.diagnostics[:timings_seconds][k] for run in runs]) for k in keys(value.diagnostics[:timings_seconds])) : Dict()
    return (; median=median(times),q25=quantile(times,0.25),q75=quantile(times,0.75),
        minimum=minimum(times), maximum=maximum(times), bytes=median(getproperty.(runs,:bytes)),
        samples=times, phase=value, repeats, stage_medians=stages)
end

function reference_representatives(f)
    f.id == "A" && return [[1],[2]]
    f.id == "B" && return [[1,1],[2,1],[1,3]]
    f.id == "C" && return [[0,1,2],[0,1,3]]
    f.id == "D" && return [[0,0,2],[0,1,1],[1,1,1]]
    f.id == "E" && return [[1,1,2],[1,2,3]]
    return vcat([ones(Int,f.n)], [ones(Int,f.n) .+ [Int(j==i) for j in 1:f.n] for i in 1:f.n])
end

function matching(result,f)
    isempty(result.atoms) && return (;atom_error=Inf,weight_error=Inf)
    matches = [argmin([maximum(abs.(v-a)) for a in f.atoms]) for v in result.atoms]
    length(unique(matches)) == length(f.atoms) || return (;atom_error=Inf,weight_error=Inf)
    return (;atom_error=maximum(maximum(abs.(v-f.atoms[i])) for (v,i) in zip(result.atoms,matches)),
        weight_error=maximum(abs(w-f.weights[i]) for (w,i) in zip(result.weights,matches)))
end

function benchmark(f)
    println("Running ",f.id," n=",f.n); flush(stdout)
    generation = measure(() -> invariant_moments(f.atoms,f.weights,f.degrees,f.t);repeats=7)
    y = generation.phase
    yf = Dict(a=>Float64(v) for (a,v) in y)
    exact = measure(() -> extract_invariant_atoms(y,f.degrees,f.t;arithmetic=:exact,basis=f.B,seed=SEED);repeats=7)
    numeric = measure(() -> extract_invariant_atoms(yf,f.degrees,f.t;seed=SEED))
    er,nr = exact.phase,numeric.phase
    @assert er.status == nr.status == :success
    @assert er.diagnostics[:exact_full_moments_consistent]
    @assert nr.rank == length(f.atoms)
    @assert nr.diagnostics[:full_representation]
    lift = measure(() -> lift_invariant_atoms(nr;symmetry=f.symmetry);repeats=7)
    orbits = lift.phase
    @assert all(o -> o.status == :lifted && o.points === nothing,orbits)
    enumerated = f.id != "F" || f.n <= 8 ? measure(() -> lift_invariant_atoms(nr;symmetry=f.symmetry,enumerate=true);repeats=3) : nothing
    checked = enumerated === nothing ? [o.representative for o in orbits] : vcat([o.points for o in enumerated.phase]...)
    obj_error = maximum(abs(objective(f,x)) for x in checked)
    match = matching(nr,f)
    @assert match.atom_error < 1e-6 && match.weight_error < 1e-6
    @assert obj_error < 1e-12
    references=reference_representatives(f)
    @assert all(x -> iszero(objective(f,x)),references)
    return Dict("id"=>f.id,"n"=>f.n,"degrees"=>f.degrees,"order"=>f.t,
        "rank"=>nr.rank,"dimension_C"=>length(nr.diagnostics[:C]),"dimension_D"=>length(nr.diagnostics[:D]),
        "dense_dimension"=>binomial(f.n+f.t,f.t),"available_moments"=>length(y),
        "required_moments"=>length(nr.diagnostics[:required_indices]),"additional_moments"=>length(nr.diagnostics[:additional_indices]),
        "original_points"=>Int(sum(o.size for o in orbits)),"orbit_sizes"=>Int[o.size for o in orbits],
        "stabilizers"=>Int[o.stabilizer for o in orbits],"representatives"=>[o.representative for o in orbits],
        "fixture_atoms"=>f.atoms,"fixture_weights"=>string.(f.weights),
        "reference_representatives"=>references,"exact_reference_objectives_zero"=>true,
        "numeric_atoms"=>nr.atoms,"numeric_weights"=>nr.weights,
        "exact_atoms"=>er.diagnostics[:exact_atoms],"exact_weights"=>er.diagnostics[:exact_weights],
        "exact_H"=>er.H,"exact_T"=>er.multiplication,"exact_basis"=>er.basis,
        "exact_full_match"=>er.diagnostics[:exact_full_moments_consistent],"numeric_full_match"=>nr.diagnostics[:full_representation],
        "rank_C"=>nr.diagnostics[:rank_C],"rank_D"=>nr.diagnostics[:rank_D],
        "atom_error"=>match.atom_error,"weight_error"=>match.weight_error,"objective_error"=>obj_error,
        "objective_checked_points"=>length(checked),"required_residual"=>nr.diagnostics[:required_residual],
        "additional_residual"=>nr.diagnostics[:additional_residual],"commutator"=>nr.diagnostics[:commutator_scaled],
        "eigenvector_residual"=>nr.diagnostics[:joint_eigenvector_residual],"eigengap"=>nr.diagnostics[:spectral_gap],
        "weight_sum_error"=>nr.diagnostics[:weight_sum_residual],"minimum_weight"=>nr.diagnostics[:minimum_weight],
        "generation"=>generation,"exact"=>exact,"numeric"=>numeric,"lift"=>lift,
        "enumeration"=>enumerated,"phase_seconds"=>nr.diagnostics[:timings_seconds])
end

function dense_matrix(f,d)
    basis = get_basis(collect(1:f.n),d)
    exponents = [ntuple(i -> count(==(i),mon),f.n) for mon in basis]
    cache = Dict{Tuple,Float64}()
    M = zeros(length(basis),length(basis))
    for i in eachindex(basis),j in i:length(basis)
        a = exponents[i] .+ exponents[j]
        value = get!(cache,a) do
            any(isodd,a) ? 0.0 : Float64(sum(w*prod(v[k]^(a[k]÷2) for k in 1:f.n) for (v,w) in zip(f.atoms,f.weights)))
        end
        M[i,j] = M[j,i] = value
    end
    return M
end

function dense_comparison()
    f=fixture("F";n=8)
    rows=Dict[]
    for d in (3,4)
        assembly=measure(() -> dense_matrix(f,d);repeats=5)
        M=assembly.phase
        diagonalization=measure(() -> eigvals(Symmetric(M));repeats=5)
        values=diagonalization.phase
        rank=count(>(sqrt(eps(Float64))*maximum(abs,values)),values)
        @assert rank == (d == 3 ? 165 : 459)
        push!(rows,Dict("order"=>d,"dimension"=>size(M,1),"rank"=>rank,
            "assembly"=>assembly,"eigensolve"=>diagonalization,"matrix_bytes"=>sizeof(M)))
    end
    return rows
end

function stress()
    rows=Dict[]
    for T in (Float64,BigFloat), exponent in (3,6,9)
        epsilon=T(10)^(-exponent)
        atoms=[[one(T),one(T)],[1+epsilon,one(T)]]
        y=invariant_moments(atoms,[T(0.5),T(0.5)],[2,2],4)
        r=extract_invariant_atoms(y,[2,2],4;seed=SEED)
        err=r.rank == 2 ? maximum(minimum(maximum(abs.(a-v)) for v in r.atoms) for a in atoms) : missing
        push!(rows,Dict("type"=>string(T),"epsilon"=>Float64(epsilon),"rank"=>r.rank,"status"=>string(r.status),
            "atom_error"=>err,"moment_residual"=>get(r.diagnostics,:required_residual,missing),
            "rank_warning"=>r.diagnostics[:rank_resolution_warning],"spectra"=>r.diagnostics[:spectra],
            "rank_thresholds"=>r.diagnostics[:rank_thresholds],"det_H_expected"=>epsilon^2/4))
    end
    return rows
end

function perturbations()
    rows=Dict[]
    for id in ("D","E"), delta in (1e-12,1e-10,1e-8,1e-6)
        f=fixture(id); y=invariant_moments(f.atoms,f.weights,f.degrees,f.t)
        counts=Dict{String,Int}(); right_rank=0; worst=0.0
        for seed in 0:9
            rng=MersenneTwister(seed)
            noisy=Dict(a => all(iszero,a) ? 1.0 : Float64(y[a])+delta*(1+abs(Float64(y[a])))*randn(rng) for a in sort!(collect(keys(y))))
            r=extract_invariant_atoms(noisy,f.degrees,f.t;seed=SEED)
            key=string(r.status); counts[key]=get(counts,key,0)+1
            r.rank == length(f.atoms) && (right_rank+=1)
            r.status == :success && (worst=max(worst,Float64(r.diagnostics[:required_residual])))
        end
        push!(rows,Dict("id"=>id,"delta"=>delta,"trials"=>10,"statuses"=>counts,"correct_rank_trials"=>right_rank,"worst_success_residual"=>worst))
    end
    return rows
end

function clean(x)
    x === nothing && return "not_run"
    x === missing && return "unresolved"
    x isa NamedTuple && return Dict(string(k)=>clean(v) for (k,v) in pairs(x) if k != :phase)
    x isa AbstractDict && return Dict(string(k)=>clean(v) for (k,v) in x)
    x isa AbstractMatrix && return [clean(collect(row)) for row in eachrow(x)]
    x isa Union{Tuple,AbstractVector} && return [clean(v) for v in x]
    x isa Rational && return string(x)
    x isa BigFloat && return Float64(x)
    x isa BigInt && return Int(x)
    x isa Symbol && return string(x)
    return x
end

fixtures=[fixture(id) for id in ("A","B","C","D","E")]
append!(fixtures,[fixture("F";n) for n in (3,5,8,10)])
cases=[benchmark(f) for f in fixtures]
negative=extract_invariant_atoms(Dict((0,)=>1//1,(1,)=>1//1,(2,)=>1//1,(3,)=>2//1),[2],3)
@assert negative.status == :inconsistent_full_moments
results=Dict("metadata"=>Dict("run_date"=>"2026-11-02","timestamp_utc"=>string(now(UTC)),
    "julia"=>string(VERSION),"cpu"=>Sys.cpu_info()[1].model,"cpu_name"=>Sys.CPU_NAME,
    "total_memory_bytes"=>Int(Sys.total_memory()),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
    "seed"=>SEED,"float_rank_rtol"=>sqrt(eps(Float64)),"float_validation_tol"=>100sqrt(eps(Float64)),
    "bigfloat_bits"=>precision(BigFloat),"runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "git_head"=>strip(read(`git rev-parse HEAD`,String)),"source_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","invariant_extraction.jl"))))),
    "cases"=>cases,"dense_comparison"=>dense_comparison(),"stress"=>stress(),"perturbations"=>perturbations(),
    "negative_control"=>Dict("status"=>string(negative.status),"atoms"=>negative.atoms,"weights"=>negative.weights,
        "required_residual"=>negative.diagnostics[:required_residual],"additional_residual"=>negative.diagnostics[:additional_residual],
        "full_representation"=>negative.diagnostics[:full_representation]))
results["metadata"]["process_peak_rss_bytes"]=Sys.maxrss()
open(joinpath(OUT,"results.toml"),"w") do io
    TOML.print(io,clean(results))
end
println("Results saved to ",joinpath(OUT,"results.toml"))
end
