# julia --project=. example/invariant_extraction_benchmarks.jl [output.csv]
# Allocations are measured by @timed; they are not a peak resident-memory metric.
using TSSOS, Printf

function benchmark_family(n)
    atoms = [ones(Int,n)]
    append!(atoms,[ones(Int,n)+3 .* [j == i for j in 1:n] for i in 1:n])
    weights = [1//2;fill(1//(2n),n)]
    degrees = fill(2,n)
    generation = @timed invariant_moments(atoms,weights,degrees,4)
    extraction = @timed extract_invariant_atoms(generation.value,degrees,4)
    result = extraction.value
    result.status == :success || error("n=$n extraction failed: $(result.status)")
    lifting = @timed lift_invariant_atoms(result)
    all(o -> o.points === nothing,lifting.value) || error("unexpected orbit enumeration")
    return (; n, invariant_atoms=result.rank,
        original_points=sum(o.size for o in lifting.value),
        invariant_dimension=length(weighted_indices(degrees,4)), dense_dimension=binomial(n+4,4),
        generation_seconds=generation.time, extraction_seconds=extraction.time,
        basis_seconds=result.diagnostics[:timings_seconds].basis_selection,
        matrix_seconds=result.diagnostics[:timings_seconds].moment_assembly,
        rank_seconds=result.diagnostics[:timings_seconds].rank_checks,
        joint_seconds=result.diagnostics[:timings_seconds].extraction,
        validation_seconds=result.diagnostics[:timings_seconds].validation,
        lifting_seconds=lifting.time, generation_bytes=generation.bytes,
        extraction_bytes=extraction.bytes, lifting_bytes=lifting.bytes,
        moment_error=Float64(result.diagnostics[:required_residual]))
end

# Tuple keys specialize on generator count; warm up every measured dimension.
for n in (3,5,8,10)
    benchmark_family(n)
end
rows = [benchmark_family(n) for n in (3,5,8,10)]
println("Family F (orbit sizes only; no original points materialized)")
println("n   atoms   original points  invariant dim  dense dim   extraction (s)  allocated (bytes)")
for row in rows
    @printf "%2d  %5d  %15d  %13d  %9d  %15.6f  %17d\n" row.n row.invariant_atoms row.original_points row.invariant_dimension row.dense_dimension row.extraction_seconds row.extraction_bytes
end
println("\nNear-collision diagnostics (rank, status, required residual)")
for T in (Float64,BigFloat), exponent in (3,6,9)
    epsilon = T(10)^(-exponent)
    y = invariant_moments([[one(T),one(T)],[1+epsilon,one(T)]],[T(0.5),T(0.5)],[2,2],4)
    result = extract_invariant_atoms(y,[2,2],4)
    println((; type=T, epsilon, rank=result.rank, status=result.status,
        error=get(result.diagnostics,:required_residual,missing), gaps=result.diagnostics[:rank_gaps]))
end
if !isempty(ARGS)
    open(ARGS[1],"w") do io
        println(io,join(string.(keys(first(rows))),','))
        for row in rows
            println(io,join(values(row),','))
        end
    end
end
