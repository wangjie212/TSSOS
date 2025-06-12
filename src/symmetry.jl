# import StarAlgebras as SA

struct VariablePermutation{V} <: SymbolicWedderburn.ByPermutations
    variables::V
end

abstract type OnMonomials <: SymbolicWedderburn.ByLinearTransformation end

function SymbolicWedderburn.action(
    a::VariablePermutation,
    g::AbstractPermutations.AbstractPermutation,
    m::AbstractMonomial,
)
    v = a.variables
    return m(v => SymbolicWedderburn.action(a, g, v))
end

# function SymbolicWedderburn.decompose(
#     p::MultivariateBases.Polynomial,
#     hom::SymbolicWedderburn.InducedActionHomomorphism,
# )
#     return [hom[p]], [1]
# end

# function SymbolicWedderburn.decompose(
#     p::SA.AlgebraElement,
#     hom::SymbolicWedderburn.InducedActionHomomorphism,
# )
#     return [hom[k] for k in SA.supp(p)],
#     [v for (_, v) in SA.nonzero_pairs(SA.coeffs(p))]
# end

function SymbolicWedderburn.decompose(
    k::DP.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials
    I = SymbolicWedderburn._int_type(hom)
    indcs = I[hom[mono] for mono in DP.monomials(k)]
    coeffs = DP.coefficients(k)

    return indcs, coeffs
end

"""
    opt,data = tssos_symmetry(pop, x, d, group; numeq=0, SymmetricConstraint=true, QUIET=false)

Compute the symmetry adapted moment-SOS relaxation for polynomial optimization problems.

# Input arguments
- `pop`: polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function tssos_symmetry(pop, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, QUIET=false, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    return tssos_symmetry_first(pop, x, d, group, numeq=numeq, action=action, semisimple=semisimple, SymmetricConstraint=SymmetricConstraint, TS=false, QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

"""
    opt,data = tssos_symmetry_first(pop, x, d, group; numeq=0, SymmetricConstraint=true, TS="block", QUIET=false)

Compute the symmetry adapted moment-SOS relaxation for polynomial optimization problems.

# Input arguments
- `pop`: polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function tssos_symmetry_first(pop, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, TS="block", QUIET=false, merge=false, md=3, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    m = length(pop) - 1 - numeq
    time = @elapsed begin
    if action === nothing
        action = VariablePermutation(x)
    end
    monos_2d = MP.monomials(x, 0:2d)
    monos_d = MP.monomials(x, 0:d)
    wedderburn = WedderburnDecomposition(Float64, group, action, monos_2d, monos_d, semisimple=semisimple)
    basis = Vector{Vector{Vector{Poly{Float64}}}}(undef, m+1)
    basis[1] = [[item'*monos_d for item in eachrow(ele.basis)] for ele in wedderburn.Uπs] # This is the symmetry adapted basis
    for i = 2:m+1
        basis[i] = Vector{Poly{Float64}}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                ind = maxdegree.(bas) .<= d - ceil(Int, maxdegree(pop[i])/2)
                if any(ind)
                    push!(basis[i], bas[ind])
                end
            end
        else
            ind = maxdegree.(monos_d) .<= d - ceil(Int, maxdegree(pop[i])/2)
            push!(basis[i], monos_d[ind])
        end
    end
    ebasis = Vector{Vector{Poly{Float64}}}(undef, numeq)
    if numeq > 0
        if SymmetricConstraint == true
            invbas = [minimum(monos_2d[item.nzind]) for item in wedderburn.invariants] # This is the basis for the invariant ring
            for i = 1:numeq
                ebasis[i] = invbas[maxdegree.(invbas) .<= 2d-maxdegree(pop[length(pop)-numeq+i])]
            end
        else
            for i = 1:numeq
                ebasis[i] = monos_2d[maxdegree.(monos_2d) .<= 2d-maxdegree(pop[length(pop)-numeq+i])]
            end
        end
    end
    tsupp = nothing
    if TS != false
        tsupp = vcat([poly_norm(p, group, action)[1] for p in pop]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_blocks(m, numeq, tsupp, pop, basis, ebasis, group, action, TS=TS, merge=merge, md=md)
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize[1]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        println("Assembling the SDP...")
    end
    optimum,tsupp,GramMat,multiplier,SDP_status = solvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action, numeq=numeq, QUIET=QUIET, 
    solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
    data = poly_basis(pop, numeq, Float64, group, action, basis, ebasis, tsupp, blocksize, blocks, eblocks, GramMat, multiplier, SDP_status)
    return optimum,data
end

"""
    opt,data = complex_tssos_symmetry(pop, x, d, group; numeq=0, SymmetricConstraint=true, ConjugateBasis=false, QUIET=false)

Compute the symmetry adapted moment-HSOS relaxation for complex polynomial optimization problems.

# Input arguments
- `pop`: complex polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function complex_tssos_symmetry(pop, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, ConjugateBasis=false, QUIET=false, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    return complex_tssos_symmetry_first(pop, x, d, group, numeq=numeq, action=action, semisimple=semisimple, SymmetricConstraint=SymmetricConstraint, ConjugateBasis=ConjugateBasis, TS=false, QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

"""
    opt,data = complex_tssos_symmetry_first(pop, x, d, group; numeq=0, SymmetricConstraint=true, ConjugateBasis=false, TS="block", QUIET=false)

Compute the symmetry adapted moment-HSOS relaxation for complex polynomial optimization problems.

# Input arguments
- `pop`: complex polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function complex_tssos_symmetry_first(pop::Vector{Poly{T}}, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, ConjugateBasis=false, TS="block", QUIET=false, merge=false, md=3, 
    dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para()) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    coe_type = T <: Real ? Float64 : ComplexF64
    m = length(pop) - 1 - numeq
    time = @elapsed begin
    if action === nothing
        action = VariablePermutation(x)
    end
    monos_d = MP.monomials(x, 0:d)
    if ConjugateBasis == true
        monos_d = filter(item -> maxdegree(item) <= d, vec([item1*conj(item2) for item1 in monos_d, item2 in monos_d]))
    end
    monos_2d = vec([item1*conj(item2) for item1 in monos_d, item2 in monos_d])
    if ConjugateBasis == true
        unique!(monos_2d)
    end
    wedderburn = WedderburnDecomposition(coe_type, group, action, monos_2d, monos_d, semisimple=semisimple)
    basis = Vector{Vector{Vector{Poly{coe_type}}}}(undef, m+1)
    basis[1] = [[item'*monos_d for item in eachrow(ele.basis)] for ele in wedderburn.Uπs] # This is the symmetry adapted basis
    for i = 2:m+1
        basis[i] = Vector{Poly{coe_type}}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                if ConjugateBasis == true
                    ind = maxdegree.(bas) .<= d - maxhalfdegree(pop[i])
                else
                    ind = maxdegree_complex.(bas) .<= d - maxdegree_complex(pop[i])
                end
                if any(ind)
                    push!(basis[i], bas[ind])
                end
            end
        else
            if ConjugateBasis == true
                ind = maxdegree.(monos_d) .<= d - maxhalfdegree(pop[i])
            else
                ind = maxdegree.(monos_d) .<= d - maxdegree_complex(pop[i])
            end
            push!(basis[i], monos_d[ind])
        end
    end
    ebasis = Vector{Vector{Poly{Float64}}}(undef, numeq)
    if numeq > 0
        if SymmetricConstraint == true
            invbas = [minimum(monos_2d[item.nzind]) for item in wedderburn.invariants] # This is the basis for the invariant ring
            for i = 1:numeq
                if ConjugateBasis == true
                    ebasis[i] = invbas[maxdegree.(invbas) .<= 2d-maxdegree(pop[length(pop)-numeq+i])]
                else
                    ebasis[i] = invbas[maxdegree_complex.(invbas) .<= d-maxdegree_complex(pop[length(pop)-numeq+i])]
                end
            end
        else
            for i = 1:numeq
                if ConjugateBasis == true
                    ebasis[i] = monos_2d[maxdegree.(monos_2d) .<= 2d-maxdegree(pop[length(pop)-numeq+i])]
                else
                    ebasis[i] = monos_2d[maxdegree_complex.(monos_2d) .<= d-maxdegree_complex(pop[length(pop)-numeq+i])]
                end
            end
        end
    end
    tsupp = nothing
    if TS != false
        tsupp = vcat([poly_norm(p, group, action)[1] for p in pop]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_blocks(m, numeq, tsupp, pop, basis, ebasis, group, action, TS=TS, merge=merge, md=md, field="complex")
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize[1]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        println("Assembling the SDP...")
    end
    optimum,tsupp,GramMat,multiplier,SDP_status = csolvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action, numeq=numeq, QUIET=QUIET, 
    coe_type=coe_type, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
    data = poly_basis(pop, numeq, coe_type, group, action, basis, ebasis, tsupp, blocksize, blocks, eblocks, GramMat, multiplier, SDP_status)
    return optimum,data
end

function tssos_symmetry_higher!(data::poly_basis; TS="block", merge=false, md=3, field="real", QUIET=false, solver="Mosek", dualize=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    pop = data.pop
    numeq = data.numeq
    group = data.group
    action = data.action
    basis = data.basis
    ebasis = data.ebasis
    ksupp = data.ksupp
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(length(pop)-1-numeq, numeq, ksupp, pop, basis, ebasis, group, action, TS=TS, merge=merge, md=md, field=field)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the TSSOS hierarchy!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize[1]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        if field == "real"
            opt,ksupp,GramMat,multiplier,SDP_status = solvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action, numeq=numeq, 
        QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
        else
            opt,ksupp,GramMat,multiplier,SDP_status = csolvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action, numeq=numeq, 
        QUIET=QUIET, coe_type=data.coe_type, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
        end
        data.blocks = blocks
        data.eblocks = eblocks
        data.blocksize = blocksize
        data.ksupp = ksupp
        data.GramMat = GramMat
        data.multiplier = multiplier
        data.SDP_status = SDP_status
    end
    return opt,data
end

function complex_tssos_symmetry_higher!(data::poly_basis; TS="block", merge=false, md=3, QUIET=false, solver="Mosek", dualize=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    return tssos_symmetry_higher!(data, TS=TS, merge=merge, md=md, field="complex", QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

function solvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action; numeq=0, QUIET=false, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    m = length(pop) - 1 - numeq
    if solver == "Mosek"
        if dualize == false
            model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
        else
            model = Model(dual_optimizer(Mosek.Optimizer))
        end
    elseif solver == "COSMO"
        model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter, "time_limit" => cosmo_setting.time_limit))
    elseif solver == "SDPT3"
        model = Model(optimizer_with_attributes(SDPT3.Optimizer))
    elseif solver == "SDPNAL"
        model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
    else
        @error "The solver is currently not supported!"
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    lambda = @variable(model)
    @objective(model, Max, lambda)
    poly = pop[1] - lambda
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, 1+m)
    pos[1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[1]))
    for (i, bas) in enumerate(basis[1])
        pos[1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1][i])
        for (l, bs) in enumerate(blocksize[1][i])
            if bs == 1
               pos[1][i][l] = @variable(model, lower_bound=0)
               poly -= bas[blocks[1][i][l][1]]^2 * pos[1][i][l]
            else
               pos[1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
               poly -= bas[blocks[1][i][l]]' * pos[1][i][l] * bas[blocks[1][i][l]]
            end
        end
    end
    for k = 1:m
        pos[k+1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k+1]))
        for (i, bas) in enumerate(basis[k+1])
            pos[k+1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1][i])
            for (l, bs) in enumerate(blocksize[k+1][i])
                if bs == 1
                    pos[k+1][i][l] = @variable(model, lower_bound=0)
                    poly -= bas[blocks[k+1][i][l][1]]^2 * pos[k+1][i][l] * pop[k+1]
                else
                    pos[k+1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    poly -= bas[blocks[k+1][i][l]]' * pos[k+1][i][l] * bas[blocks[k+1][i][l]] * pop[k+1]
                end
            end
        end
    end
    if numeq > 0
        mul = Vector{Vector{VariableRef}}(undef, numeq)
        for k = 1:numeq
            mul[k] = @variable(model, [1:length(eblocks[k])])
            poly -= ebasis[k][eblocks[k]]' * mul[k] * pop[k+m+1]
        end
    end
    tsupp,coe = poly_norm(poly, group, action)
    @constraint(model, coe .== 0)
    if QUIET == false
        println("Solving the SDP...")
    end
    time = @elapsed begin
    optimize!(model)
    end
    if QUIET == false
        println("SDP solving time: $time seconds.")
    end
    SDP_status = termination_status(model)
    if SDP_status != MOI.OPTIMAL
        println("termination status: $SDP_status")
        status = primal_status(model)
        println("solution status: $status")
    end
    optimum = objective_value(model)
    @show optimum
    GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, m+1)
    multiplier = nothing
    for i = 1:m+1
        GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, length(basis[i]))
        for j = 1:length(basis[i])
            GramMat[i][j] = [value.(pos[i][j][k]) for k = 1:cl[i][j]]
        end
        if numeq > 0
            multiplier = [value.(mul[j]) for j = 1:numeq]
        end
    end
    return optimum,tsupp,GramMat,multiplier,SDP_status
end

function csolvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action; numeq=0, QUIET=false, coe_type=ComplexF64, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    m = length(pop) - 1 - numeq
    if solver == "Mosek"
        if dualize == false
            model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
        else
            model = Model(dual_optimizer(Mosek.Optimizer))
        end
    elseif solver == "COSMO"
        model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter, "time_limit" => cosmo_setting.time_limit))
    elseif solver == "SDPT3"
        model = Model(optimizer_with_attributes(SDPT3.Optimizer))
    elseif solver == "SDPNAL"
        model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
    else
        @error "The solver is currently not supported!"
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    lambda = @variable(model)
    @objective(model, Max, lambda)
    poly = pop[1] - lambda
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, 1+m)
    pos[1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[1]))
    for (i, bas) in enumerate(basis[1])
        pos[1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1][i])
        for (l, bs) in enumerate(blocksize[1][i])
            if bs == 1
                pos[1][i][l] = @variable(model, lower_bound=0)
                poly -= bas[blocks[1][i][l][1]] * pos[1][i][l] * conj(bas[blocks[1][i][l][1]])
            else
                if coe_type == ComplexF64
                    pos[1][i][l] = @variable(model, [1:2bs, 1:2bs], PSD)
                    poly -= bas[blocks[1][i][l]]' * (pos[1][i][l][1:bs,1:bs]+pos[1][i][l][bs+1:2bs,bs+1:2bs]+pos[1][i][l][bs+1:2bs,1:bs]*im-pos[1][i][l][1:bs,bs+1:2bs]*im) * bas[blocks[1][i][l]]
                else
                    pos[1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    poly -= bas[blocks[1][i][l]]' * pos[1][i][l] * bas[blocks[1][i][l]]
                end
            end
        end
    end
    for k = 1:m
        pos[k+1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k+1]))
        for (i, bas) in enumerate(basis[k+1])
            pos[k+1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1][i])
            for (l, bs) in enumerate(blocksize[k+1][i])
                if bs == 1
                    pos[k+1][i][l] = @variable(model, lower_bound=0)
                    poly -= bas[blocks[k+1][i][l][1]] * conj(bas[blocks[k+1][i][l][1]]) * pos[k+1][i][l] * pop[k+1]
                else
                    if coe_type == ComplexF64
                        pos[k+1][i][l] = @variable(model, [1:2bs, 1:2bs], PSD)
                        poly -= bas[blocks[k+1][i][l]]' * (pos[k+1][i][l][1:bs,1:bs]+pos[k+1][i][l][bs+1:2bs,bs+1:2bs]+pos[k+1][i][l][bs+1:2bs,1:bs]*im-pos[k+1][i][l][1:bs,bs+1:2bs]*im) * bas[blocks[k+1][i][l]] * pop[k+1]
                    else
                        pos[k+1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                        poly -= bas[blocks[k+1][i][l]]' * pos[k+1][i][l] * bas[blocks[k+1][i][l]] * pop[k+1]
                    end
                end
            end
        end
    end
    if numeq > 0
        mul = Vector{Vector{VariableRef}}(undef, numeq)
        for k = 1:numeq
            temp = ebasis[k][eblocks[k]][[item <= conj(item) for item in ebasis[k][eblocks[k]]]]
            lb = length(temp)
            if coe_type == ComplexF64
                mul[k] = @variable(model, [1:2*lb])
                tau = sum(temp .* (mul[k][1:lb] .+ mul[k][lb+1:2*lb]*im))
            else
                mul[k] = @variable(model, [1:lb])
                tau = sum(temp .* mul[k])
            end
            poly -= (tau + conj(tau)) * pop[k+m+1]
        end
    end
    tsupp,coe = poly_norm(poly, group, action)
    coe = complex_reduce(tsupp, coe)[2]
    @constraint(model, coe .== 0)
    if QUIET == false
        println("Solving the SDP...")
    end
    time = @elapsed begin
    optimize!(model)
    end
    if QUIET == false
        println("SDP solving time: $time seconds.")
    end
    SDP_status = termination_status(model)
    if SDP_status != MOI.OPTIMAL
        println("termination status: $SDP_status")
        status = primal_status(model)
        println("solution status: $status")
    end
    optimum = objective_value(model)
    @show optimum
    GramMat = Vector{Vector{Vector{Union{Float64,Matrix{coe_type}}}}}(undef, m+1)
    multiplier = nothing
    for i = 1:m+1
        GramMat[i] = Vector{Vector{Union{Float64,Matrix{coe_type}}}}(undef, length(basis[i]))
        for j = 1:length(basis[i])
            if coe_type == Float64
                GramMat[i][j] = [value.(pos[i][j][k]) for k = 1:cl[i][j]]
            else
                GramMat[i][j] = Vector{Union{Float64,Matrix{coe_type}}}(undef, cl[i][j])
                for k = 1:cl[i][j]
                    bs = blocksize[i][j][k]
                    if bs == 1
                        GramMat[i][j][k] = value(pos[i][j][k])
                    else
                        GramMat[i][j][k] = value.(pos[i][j][k][1:bs,1:bs])+value.(pos[i][j][k][bs+1:2*bs,bs+1:2*bs])+value.(pos[i][j][k][bs+1:2*bs,1:bs])*im-value.(pos[i][j][k][1:bs,bs+1:2*bs])*im
                    end
                end
            end
        end
    end
    if numeq > 0
        multiplier = Vector{Poly{coe_type}}(undef, numeq)
        for k = 1:numeq
            temp = ebasis[k][eblocks[k]][[item <= conj(item) for item in ebasis[k][eblocks[k]]]]
            lb = length(temp)
            if coe_type == Float64
                tau = sum(temp .* value.(mul[k]))
            else
                tau = sum(temp .* (value.(mul[k][1:lb]) .+ value.(mul[k][lb+1:2*lb])*im))
            end
            multiplier[k] = tau + conj(tau)
        end
    end
    return optimum,tsupp,GramMat,multiplier,SDP_status
end

"""
    supp,coe = poly_norm(p, group, action)

Compute the normal form of a polynomial under the action of a group.

# Input arguments
- `p`: polynomial
- `group`: group
- `action`: describing how the group acts on monomials

# Output arguments
- `supp`: support of the normal form
- `coe`: coefficients of  the normal form
"""
function poly_norm(p, group, action)
    supp = MP.monomials(p)
    coe = MP.coefficients(p)
    nsupp = [normalform(mon, group, action) for mon in supp]
    supp = copy(nsupp)
    unique!(supp)
    sort!(supp)
    ncoe = zeros(typeof(coe[1]), length(supp))
    for i = 1:length(nsupp)
        @inbounds ind = bfind(supp, length(supp), nsupp[i])
        @inbounds ncoe[ind] += coe[i]
    end
    return supp, ncoe
end

function complex_reduce(supp, coe)
    nsupp = [supp[1]]
    ncoe = [coe[1]]
    for i = 2:length(supp)
        if conj(supp[i]) == supp[i] || (bfind(nsupp, length(nsupp), supp[i]) === nothing && bfind(nsupp, length(nsupp), conj(supp[i])) === nothing)
            push!(nsupp, supp[i])
            push!(ncoe, coe[i])
        end
    end
    return nsupp, ncoe
end

"""
    mon = normalform(mon, group, action)

Compute the normal form of a monomial under the action of a group.

# Input arguments
- `mon`: monomial
- `group`: group
- `action`: describing how the group acts on monomials

# Output arguments
- `mon`: normal form of the monomial
"""
function normalform(mon, group, action)
    return minimum([SymbolicWedderburn.action(action, g, mon) for g in group])
end

function get_blocks(m::Int, l::Int, tsupp, pop, basis::Vector{Vector{Vector{Poly{T}}}}, ebasis::Vector{Vector{Poly{Float64}}}, 
    group, action; TS="block", merge=false, md=3, field="real") where {T<:Number}
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    blocksize = Vector{Vector{Vector{Int}}}(undef, m+1)
    cl = Vector{Vector{Int}}(undef, m+1)
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [[Vector(1:length(basis[k][i]))] for i = 1:length(basis[k])],[[length(basis[k][i])] for i = 1:length(basis[k])],ones(Int, length(basis[k]))
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else
        for k = 1:m+1
            blocks[k] = Vector{Vector{Vector{Int}}}(undef, length(basis[k]))
            blocksize[k] = Vector{Vector{Int}}(undef, length(basis[k]))
            cl[k] = Vector{Int}(undef, length(basis[k]))
            for (i, ba) in enumerate(basis[k])
                if k == 1
                    G = get_graph(tsupp, ba, group, action, field=field)
                else
                    G = get_graph(tsupp, ba, group, action, g=pop[k-1], field=field)
                end
                if TS == "block"
                    blocks[k][i] = connected_components(G)
                    blocksize[k][i] = length.(blocks[k][i])
                    cl[k][i] = length(blocksize[k][i])
                else
                    blocks[k][i],cl[k][i],blocksize[k][i] = chordal_cliques!(G, method=TS)
                    if merge == true
                        blocks[k][i],cl[k][i],blocksize[k][i] = clique_merge!(blocks[k], d=md, QUIET=true)
                    end
                end
            end
        end
        for k = 1:l
            eblocks[k] = get_eblock(tsupp, pop[k+m+1], ebasis[k], group, action)
        end
    end
    return blocks,cl,blocksize,eblocks
end

function get_graph(tsupp, basis::Vector{Poly{T}}, group, action; g=1, field="real") where {T<:Number}
    lb = length(basis)
    G = SimpleGraph(lb)
    ltsupp = length(tsupp)
    for i = 1:lb, j = i+1:lb
        flag = 0
        for item in poly_norm(field == "real" ? basis[i] * basis[j] * g : basis[i] * conj(basis[j]) * g, group, action)[1]
            if lbfind(tsupp, ltsupp, item) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_eblock(tsupp, h, basis::Vector{Poly{Float64}}, group, action)
    ltsupp = length(tsupp)
    eblock = Int[]
    for (i, ba) in enumerate(basis)
        flag = 0
        for item in poly_norm(ba * h, group, action)[1]
            if lbfind(tsupp, ltsupp, item) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            push!(eblock, i)
        end
    end
    return eblock
end

function add_psatz_symmetry!(model, nonneg::DP.Polynomial{V, M, T}, vars, ineq_cons, eq_cons, order, group; action=nothing, semisimple=false, SymmetricConstraint=true, TS="block", SO=1, merge=false, md=3, QUIET=false) where {V, M, T<:Union{Number,AffExpr}}
    m = length(ineq_cons)
    l = length(eq_cons)
    if action === nothing
        action = VariablePermutation(vars)
    end
    monos_2d = MP.monomials(vars, 0:2*order)
    monos_d = MP.monomials(vars, 0:order)
    wedderburn = WedderburnDecomposition(Float64, group, action, monos_2d, monos_d, semisimple=semisimple)
    basis = Vector{Vector{Vector{Poly{Float64}}}}(undef, m+1)
    basis[1] = [[item'*monos_d for item in eachrow(ele.basis)] for ele in wedderburn.Uπs]
    for i = 1:m
        basis[i+1] = Vector{Poly{Float64}}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                ind = maxdegree.(bas) .<= order - ceil(Int, maxdegree(ineq_cons[i])/2)
                if any(ind)
                    push!(basis[i+1], bas[ind])
                end
            end
        else
            ind = maxdegree.(monos_d) .<= order - ceil(Int, maxdegree(ineq_cons[i])/2)
            push!(basis[i+1], monos_d[ind])
        end
    end
    ebasis = Vector{Vector{Poly{Float64}}}(undef, l)
    if l > 0
        if SymmetricConstraint == true
            invbas = [minimum(monos_2d[item.nzind]) for item in wedderburn.invariants] # This is the basis for the invariant ring
            for i = 1:l
                ebasis[i] = invbas[maxdegree.(invbas) .<= 2*order-maxdegree(eq_cons[i])]
            end
        else
            for i = 1:l
                ebasis[i] = monos_2d[maxdegree.(monos_2d) .<= 2*order-maxdegree(eq_cons[i])]
            end
        end
    end
    tsupp = nothing
    if TS != false
        pop = [nonneg]
        if !isempty(ineq_cons)
            pop = [pop; ineq_cons]
        end
        if !isempty(eq_cons)
            pop = [pop; eq_cons]
        end
        tsupp = vcat([poly_norm(p, group, action)[1] for p in pop]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_blocks(m, l, tsupp, ineq_cons, eq_cons, basis, ebasis, group, action, TS=TS, SO=SO, merge=merge, md=md, QUIET=QUIET)
    poly = nonneg
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, 1+m)
    pos[1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[1]))
    for (i, bas) in enumerate(basis[1])
        pos[1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1][i])
        for (l, bs) in enumerate(blocksize[1][i])
            if bs == 1
               pos[1][i][l] = @variable(model, lower_bound=0)
               poly -= bas[blocks[1][i][l][1]]^2 * pos[1][i][l]
            else
               pos[1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
               poly -= bas[blocks[1][i][l]]' * pos[1][i][l] * bas[blocks[1][i][l]]
            end
        end
    end
    for k = 1:m
        pos[k+1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k+1]))
        for (i, bas) in enumerate(basis[k+1])
            pos[k+1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1][i])
            for (l, bs) in enumerate(blocksize[k+1][i])
                if bs == 1
                    pos[k+1][i][l] = @variable(model, lower_bound=0)
                    poly -= bas[blocks[k+1][i][l][1]]^2 * pos[k+1][i][l] * ineq_cons[k]
                else
                    pos[k+1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    poly -= bas[blocks[k+1][i][l]]' * pos[k+1][i][l] * bas[blocks[k+1][i][l]] * ineq_cons[k]
                end
            end
        end
    end
    mul = nothing
    if l > 0
        mul = Vector{Vector{VariableRef}}(undef, l)
        for k = 1:l
            mul[k] = @variable(model, [1:length(eblocks[k])])
            poly -= ebasis[k][eblocks[k]]' * mul[k] * eq_cons[k]
        end
    end
    tsupp,coe = poly_norm(poly, group, action)
    @constraint(model, coe .== 0)
    info = poly_basis(nothing, l, nothing, group, action, basis, ebasis, tsupp, blocksize, blocks, eblocks, pos, mul, nothing)
    return info
end

function get_blocks(m::Int, l::Int, tsupp, ineq_cons, eq_cons, basis, ebasis, group, action; TS="block", SO=1, merge=false, md=3, QUIET=false)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    blocksize = Vector{Vector{Vector{Int}}}(undef, m+1)
    cl = Vector{Vector{Int}}(undef, m+1)
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [[Vector(1:length(item))] for item in basis[k]],[[length(item)] for item in basis[k]],ones(Int, length(basis[k]))
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for k = 1:m+1
                blocks[k] = Vector{Vector{Vector{Int}}}(undef, length(basis[k]))
                blocksize[k] = Vector{Vector{Int}}(undef, length(basis[k]))
                cl[k] = Vector{Int}(undef, length(basis[k]))
                for (i, ba) in enumerate(basis[k])
                    if k == 1
                        G = get_graph(tsupp, ba, group, action)
                    else
                        G = get_graph(tsupp, ba, g=ineq_cons[k-1], group, action)
                    end
                    if TS == "block"
                        blocks[k][i] = connected_components(G)
                        blocksize[k][i] = length.(blocks[k][i])
                        cl[k][i] = length(blocksize[k][i])
                    else
                        blocks[k][i],cl[k][i],blocksize[k][i] = chordal_cliques!(G, method=TS)
                        if merge == true
                            blocks[k][i],cl[k][i],blocksize[k][i] = clique_merge!(blocks[k], d=md, QUIET=true)
                        end
                    end
                end
            end
            for k = 1:l
                eblocks[k] = get_eblock(tsupp, eq_cons[k], ebasis[k], group, action)
            end
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                println("No higher TS step of the TSSOS hierarchy!")
                break
            end
            if i < SO
                tsupp = DP.Monomial[]
                for t = 1:length(blocks[1]), s = 1:length(blocksize[1][t]), j = 1:blocksize[1][t][s], r = j:blocksize[1][t][s]  
                    append!(tsupp, poly_norm(basis[1][t][blocks[1][t][s][j]]*basis[1][t][blocks[1][t][s][r]], group, action)[1])
                end
                unique!(tsupp)
                sort!(tsupp)
            end
        end
    end
    return blocks,cl,blocksize,eblocks
end
