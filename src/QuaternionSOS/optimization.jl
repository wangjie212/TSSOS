"""
    qtssos(pop, z, n, d; numeq=0, CS="MF", TS=false, optimizer=nothing, kwargs...)

Compute a quaternion SOS lower bound at order `d`. The noncommuting variables
`z[1:n]` are followed by their conjugates `z[n+1:2n]`. Put the Hermitian objective
first and the `numeq` equalities last. Supply an SDP optimizer or an empty JuMP
`model`; otherwise TSSOS selects its default solver. `solve=false` returns
`nothing` after computing supports, without selecting a solver.

See `docs/src/quaternion.md` for sparsity, unit-norm, and coefficient options.
"""
function qtssos(pop, z, n::Int, d; fsupp=nothing, fcoe=nothing, numeq=0, rncnumeq=0, qncnumeq=0, nb=0, CS="MF", cliques=[], TS=false,
    merge=false, md=3, optimizer=nothing, model=nothing, QUIET=false, solve=true, ipart=true, conjubasis=true, addrcons=true, addicons=false)
    n > 0 && length(z) == 2n || throw(ArgumentError("z must contain n variables followed by their n conjugates."))
    length(unique(z)) == 2n || throw(ArgumentError("Quaternion variables must be distinct."))
    2n <= typemax(UInt16) || throw(ArgumentError("Too many quaternion variables."))
    d isa Integer && d >= 0 || throw(ArgumentError("d must be a nonnegative integer."))
    isempty(pop) && throw(ArgumentError("pop must contain an objective."))
    0 <= numeq < length(pop) || throw(ArgumentError("numeq must count trailing constraints."))
    0 <= rncnumeq && 0 <= qncnumeq && rncnumeq + qncnumeq <= numeq ||
        throw(ArgumentError("Noncommutative equality counts must sum to at most numeq."))
    0 <= nb <= n || throw(ArgumentError("nb must lie between zero and n."))
    qncnumeq > 0 && !ipart && throw(ArgumentError("Quaternion equalities require ipart=true."))
    (fsupp === nothing) == (fcoe === nothing) || throw(ArgumentError("fsupp and fcoe must be provided together."))
    if fsupp !== nothing
        length(fsupp) == length(fcoe) || throw(ArgumentError("fsupp and fcoe must have equal lengths."))
    end
    input_constraints = length(pop) - 1
    if addrcons
        rcomm_constraints = generate_rcomm_constraints(z, n)
        # Quaternion equalities remain the final qncnumeq constraints.
        split = length(pop) - qncnumeq
        supp,coe = QPolys_info(vcat(pop[1:split], rcomm_constraints, pop[split+1:end]), z, n)
        numeq = numeq + length(rcomm_constraints)
        rncnumeq = rncnumeq + length(rcomm_constraints)
    else
        supp,coe = QPolys_info(pop, z, n)
    end
    if ipart && fsupp !== nothing && fcoe !== nothing
        supp[1] = standardterm.(fsupp, n)
        coe[1] = QuaternionF64.(fcoe)
    end
    if !QUIET
        println("*********************************** TSSOS ***********************************")
        println("TSSOS is launching...")
    end
    if nb > 0
        supp[1],coe[1] = qresort(supp[1],coe[1],n;nb=nb)
    end
    m = length(supp) - 1
    dc = zeros(Int, m)
    for i = 1:m
        dc[i] = maximum([length(supp[i+1][j][1])+length(supp[i+1][j][2])+length(supp[i+1][j][3]) for j=1:length(supp[i+1])]; init=0)
    end
    maximum((sum(length, term) for support in supp for term in support); init=0) <= 2d ||
        throw(ArgumentError("The relaxation order is too small for the polynomial degrees."))
    if cliques != []
        cliques = [sort!(unique(Int.(clique))) for clique in cliques]
        all(v -> 1 <= v <= n, Iterators.flatten(cliques)) || throw(ArgumentError("Invalid variable in cliques."))
        sort!(unique(vcat(cliques...))) == collect(1:n) || throw(ArgumentError("cliques must cover all variables."))
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        time = @elapsed begin
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp(n, m, rncnumeq, supp, alg=CS, qncnumeq=qncnumeq, QUIET=QUIET)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
        end
    end
    I,J,ncc = assign_constraint(n, m, numeq, supp, cliques, cql)
    # Commutation identities between separate cliques need no local multiplier.
    any(k -> k <= input_constraints-qncnumeq || k > m-qncnumeq, ncc) && throw(ArgumentError("Every input constraint must be contained in a variable clique; enlarge cliques or use CS=false."))
    rlorder = d*ones(Int, cql)
    ebasis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    basis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    for i = 1:cql
        ebasis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(J[i]))
        basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1)
        basis[i][1] = get_qncbasis(cliques[i], n, d; conjubasis=conjubasis)
        if nb > 0
            basis[i][1] = qreduce_unitnorm.(basis[i][1], n, nb=nb)
            unique!(basis[i][1])
        end
        for s = 1:length(I[i])
            basis[i][s+1] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[I[i][s]]/2)))
            if nb > 0
                basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], n, nb=nb)
                unique!(basis[i][s+1])
            end
        end
        for s = 1:length(J[i])
            ebasis[i][s] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[J[i][s]]/2)))
            if nb > 0
                ebasis[i][s] = qreduce_unitnorm.(ebasis[i][s], n, nb=nb)
                unique!(ebasis[i][s])
            end
        end
    end
    tsupp = copy(supp[1])
    for i = 2:m+1, j = 1:length(supp[i])
        push!(tsupp, supp[i][j])
    end
    sort!(tsupp)
    unique!(tsupp)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed blocks,eblocks,cl,blocksize = get_blocks(n,rlorder, I, J, supp, cliques, cql, tsupp, basis, ebasis, TS=TS, merge=merge, md=md, nb=nb)
    if QUIET == false
        mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt, ksupp, SDP_status =  qsolvesdp(n, m, rncnumeq, qncnumeq, supp, coe, basis, ebasis, cql, I, J, ncc, blocks, eblocks, cl, blocksize, QUIET=QUIET, TS=TS, optimizer=optimizer, model=model, solve=solve, ipart=ipart,
    nb=nb, addicons=addicons)
    return opt
end

function QPolys_info(pop, z, n)
    polys = [p isa Number ? p + 0.0*z[1] : p for p in pop]
    coe = Vector{QuaternionF64}[QuaternionF64.(coefficients(p)) for p in polys]
    supp = Vector{Vector{Vector{UInt16}}}[
        Vector{Vector{UInt16}}[mono_to_term(m, z, n) for m in monomials(p)] for p in polys]
    return supp, coe
end

function qsolvesdp(n, m, rncnumeq, qncnumeq, supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, basis, ebasis, cql, I, J, ncc, blocks, eblocks, cl, blocksize;
    nb=0, QUIET=false, TS=false, optimizer=nothing, model=nothing, solve=true,  ipart=true,
    addicons=false)
    tsupp = Vector{Vector{UInt16}}[]
    ttsupp = Vector{Vector{UInt16}}[]
    if addicons
        addsupp,addcoe = generate_icomm_constraints(n)
    end
    for i = 1:cql
        a = 1
        for s = 1 : a, j = 1:cl[i][s], k = 1:blocksize[i][s][j], r = 1:blocksize[i][s][j]
            at=basis[i][s][blocks[i][s][j][k]]
            b=basis[i][s][blocks[i][s][j][r]]
            if ipart
                @inbounds bi = qtermadd(b,at,n)
            else
                @inbounds bi = qtermaddleft(at,b,n)
            end
            bi, bi_cano = reduced_terms(bi, n, nb, ipart)
            push!(tsupp, bi)
            if ipart
                push!(ttsupp,bi_cano)
            end
        end
        for (j, k) in enumerate(J[i])
            if m - rncnumeq - qncnumeq < J[i][j] <= m - qncnumeq
                bs = length(eblocks[i][j])
                for t = 1:bs, r = 1:bs
                    for s = 1:length(supp[k+1])
                        a=ebasis[i][j][eblocks[i][j][t]]
                        b=supp[k+1][s]
                        c=ebasis[i][j][eblocks[i][j][r]]
                        if ipart
                            @inbounds bi = qtermadd3(b,c,a,n)
                        else
                            @inbounds bi = qtermadd3left(a,b,c,n)
                        end
                        bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                        push!(tsupp, bi)
                        if ipart
                            push!(ttsupp,bi_cano)
                        end
                    end
                    for s = 1:length(supp[k+1])
                        a=ebasis[i][j][eblocks[i][j][t]]
                        b=supp[k+1][s]
                        c=ebasis[i][j][eblocks[i][j][r]]
                        a_star = star(a,n)
                        b_star = star(b,n)
                        if ipart
                            @inbounds bi = qtermadd3(c,a_star,b,n)
                        else
                            @inbounds bi = qtermadd3left(a,b_star,c,n)
                        end
                        bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                        push!(tsupp, bi)
                        if ipart
                            push!(ttsupp,bi_cano)
                        end
                    end
                    if addicons
                        for ss = 1:length(addsupp)
                            for tt = 1:4
                                a=ebasis[i][j][eblocks[i][j][t]]
                                b=ebasis[i][j][eblocks[i][j][r]]
                                c=addsupp[ss][tt][1]
                                d=addsupp[ss][tt][2]
                                a_star = star(a,n)
                                @inbounds bi = qtermadd4(d,b,a_star,c,n)
                                bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                push!(tsupp, bi)
                                if ipart
                                    push!(ttsupp,bi_cano)
                                end
                            end
                        end
                        for ss = 1:length(addsupp)
                            for tt = 1:4
                                a=ebasis[i][j][eblocks[i][j][t]]
                                b=ebasis[i][j][eblocks[i][j][r]]
                                c=addsupp[ss][tt][1]
                                d=addsupp[ss][tt][2]
                                a_star = star(a,n)
                                c_star = star(c,n)
                                d_star = star(d,n)
                                @inbounds bi = qtermadd4(c_star,b,a_star,d_star,n)
                                bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                push!(tsupp, bi)
                                if ipart
                                    push!(ttsupp,bi_cano)
                                end
                            end
                        end
                    end
                end
            elseif J[i][j] > m - qncnumeq
                selected = ebasis[i][j][eblocks[i][j]]
                for a in selected, b in selected, hword in supp[k+1]
                    word = qtermadd3(hword, b, a, n)
                    for term in (word, star(word, n))
                        term, cyclic = reduced_terms(term, n, nb, true)
                        push!(tsupp, term)
                        push!(ttsupp, cyclic)
                    end
                end
            end
        end
    end
    if TS != false
        gsupp = get_gsupp(n, basis, ebasis, supp, cql, I, J, ncc, blocks, eblocks, cl, blocksize, nb=nb)
        for term in gsupp
            word, cyclic = reduced_terms(term, n, nb, ipart)
            push!(tsupp, word)
            ipart && push!(ttsupp, cyclic)
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    sort!(ttsupp)
    unique!(ttsupp)
    ksupp = tsupp
    objv = SDP_status= nothing
    if solve == true
        ltsupp = length(tsupp)
        lttsupp = length(ttsupp)
        if QUIET == false
            println("Assembling the SDP...")
        end
        model = _sdp_model(model, optimizer, false)
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        rcons = [AffExpr(0) for i=1:(ipart ? lttsupp : ltsupp)]
        if ipart==true || qncnumeq > 0
            icons = [AffExpr(0) for i=1:ltsupp]
            jcons = [AffExpr(0) for i=1:ltsupp]
            kcons = [AffExpr(0) for i=1:ltsupp]
        end
        pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, cql)
        for i = 1:cql
            a = 1
            pos[i] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(I[i])+a)
            for j = 1:a
                pos[i][j] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][j])
                for l = 1:cl[i][j]
                    @inbounds bs = blocksize[i][j][l]
                    if bs == 1
                        pos[i][j][l] = @variable(model, lower_bound=0)
                        at=basis[i][j][blocks[i][j][l][1]]
                        b=basis[i][j][blocks[i][j][l][1]]
                        if ipart
                            @inbounds bi = qtermadd(b,at,n)
                        else
                            @inbounds bi = qtermaddleft(at,b,n)
                        end
                        bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                        Locb = bfind(tsupp, bi)
                        if ipart == true
                            Locb_cano = bfind(ttsupp, bi_cano)
                            @inbounds add_to_expression!(rcons[Locb_cano], pos[i][j][l])
                        else
                            @inbounds add_to_expression!(rcons[Locb], pos[i][j][l])
                        end
                    else
                        if ipart == true
                            pos[i][j][l] = @variable(model, [1:4bs, 1:4bs], PSD)
                        else
                            pos[i][j][l] = @variable(model, [1:bs, 1:bs], PSD)
                        end
                        for t = 1:bs, r = 1:bs
                            at=basis[i][j][blocks[i][j][l][t]]
                            b=basis[i][j][blocks[i][j][l][r]]
                            if ipart
                                @inbounds bi = qtermadd(b,at,n)
                            else
                                @inbounds bi = qtermaddleft(at,b,n)
                            end
                            bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                            Locb = bfind(tsupp, bi)
                            if ipart == true
                                Locb_cano = bfind(ttsupp, bi_cano)
                                @inbounds add_to_expression!(rcons[Locb_cano], pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs]+pos[i][j][l][t+2*bs,r+2*bs]+pos[i][j][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], pos[i][j][l][t+bs,r]-pos[i][j][l][t,r+bs]+pos[i][j][l][t+3*bs,r+2*bs]-pos[i][j][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(jcons[Locb], pos[i][j][l][t+2*bs,r]-pos[i][j][l][t,r+2*bs]-pos[i][j][l][t+3*bs,r+bs]+pos[i][j][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], pos[i][j][l][t+3*bs,r]-pos[i][j][l][t,r+3*bs]+pos[i][j][l][t+2*bs,r+bs]-pos[i][j][l][t+bs,r+2*bs])
                            else
                                @inbounds add_to_expression!(rcons[Locb], pos[i][j][l][t,r])
                            end
                        end
                    end
                end
            end
        end
        for i = 1:cql, (j, k) in enumerate(I[i])
            a = 1
            pos[i][j+a] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][j+a])
            for l = 1: cl[i][j+a]
                bs = blocksize[i][j+a][l]
                if bs == 1
                    pos[i][j+a][l] = @variable(model, lower_bound=0)
                    for s = 1:length(supp[k+1])
                        at=basis[i][j+a][blocks[i][j+a][l][1]]
                        b=supp[k+1][s]
                        c=basis[i][j+a][blocks[i][j+a][l][1]]
                        if ipart
                            @inbounds bi = qtermadd3(b,c,at,n)
                        else
                            @inbounds bi = qtermadd3left(at,b,c,n)
                        end
                        bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                        Locb = bfind(tsupp, bi)
                        if ipart == true
                            Locb_cano = bfind(ttsupp, bi_cano)
                            @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), pos[i][j+a][l])
                            @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], pos[i][j+a][l])
                            @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], pos[i][j+a][l])
                            @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[i][j+a][l])
                        else
                            @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[i][j+a][l])
                        end
                    end
                else
                    if ipart == true
                        pos[i][j+a][l]= @variable(model, [1:4bs, 1:4bs], PSD)
                    else
                        pos[i][j+a][l]= @variable(model, [1:bs, 1:bs], PSD)
                    end
                    for t = 1:bs, r = 1:bs
                        for s = 1:length(supp[k+1])
                            at=basis[i][j+a][blocks[i][j+a][l][t]]
                            b=supp[k+1][s]
                            c=basis[i][j+a][blocks[i][j+a][l][r]]
                            if ipart
                                @inbounds bi = qtermadd3(b,c,at,n)
                            else
                                @inbounds bi = qtermadd3left(at,b,c,n)
                            end
                            bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                            Locb = bfind(tsupp, bi)
                            if ipart == true
                                Locb_cano = bfind(ttsupp, bi_cano)
                                @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])

                                @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])

                                @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                            else
                                @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t,r])
                            end
                        end
                    end
                end
            end
        end

        epos = Vector{Vector{Any}}(undef, cql)
        if addicons
            addepos = Vector{Vector{Any}}(undef, cql)
        end
        for i = 1:cql
            if !isempty(J[i])
                epos[i] = Vector{Any}(undef, length(J[i]))
                if addicons
                    if cql > 1
                        addepos[i] = Vector{Any}(undef, 3*n*(n-1))
                    else
                        addepos[i] = Vector{Any}(undef, 3*length(J[i]))
                    end
                end
                for (j, k) in enumerate(J[i])
                    if J[i][j] <= m - rncnumeq - qncnumeq
                        bs = length(eblocks[i][j])
                        if bs == 1
                            epos[i][j] = @variable(model)
                            for s = 1:length(supp[k+1])
                                a=ebasis[i][j][eblocks[i][j][1]]
                                b=supp[k+1][s]
                                c=ebasis[i][j][eblocks[i][j][1]]
                                if ipart
                                    @inbounds bi = qtermadd3(b,c,a,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b,c,n)
                                end
                                bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                Locb = bfind(tsupp, bi)
                                if ipart == true
                                    Locb_cano = bfind(ttsupp, bi_cano)
                                    @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j])
                                end
                            end
                        else
                            if ipart == true
                                epos[i][j]= @variable(model, [1:4bs, 1:4bs], Symmetric)
                            else
                                epos[i][j]= @variable(model, [1:bs, 1:bs], Symmetric)
                            end
                            for t = 1:bs, r = 1:bs
                                for s = 1:length(supp[k+1])
                                    a=ebasis[i][j][eblocks[i][j][t]]
                                    b=supp[k+1][s]
                                    c=ebasis[i][j][eblocks[i][j][r]]
                                    if ipart
                                        @inbounds bi = qtermadd3(b,c,a,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b,c,n)
                                    end
                                    bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                    Locb = bfind(tsupp, bi)
                                    if ipart == true
                                        Locb_cano = bfind(ttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])

                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][t,r])
                                    end
                                end
                            end
                        end
                    elseif J[i][j] <= m - qncnumeq
                        bs = length(eblocks[i][j])
                        if bs == 1
                            epos[i][j] = @variable(model)
                            for s = 1:length(supp[k+1])
                                a=ebasis[i][j][eblocks[i][j][1]]
                                b=supp[k+1][s]
                                c=ebasis[i][j][eblocks[i][j][1]]
                                if ipart
                                    @inbounds bi = qtermadd3(b,c,a,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b,c,n)
                                end
                                bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                Locb = bfind(tsupp, bi)
                                if ipart == true
                                    Locb_cano = bfind(ttsupp, bi_cano)
                                    @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j])
                                end
                            end
                            for s = 1:length(supp[k+1])
                                a=ebasis[i][j][eblocks[i][j][1]]
                                b=supp[k+1][s]
                                c=ebasis[i][j][eblocks[i][j][1]]
                                a_star = star(a,n)
                                b_star = star(b,n)
                                if ipart
                                    @inbounds bi = qtermadd3(c,a_star,b,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b_star,c,n)
                                end
                                bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                Locb = bfind(tsupp, bi)
                                if ipart == true
                                    Locb_cano = bfind(ttsupp, bi_cano)
                                    @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j])
                                end
                            end
                        else
                            if ipart == true
                                epos[i][j]= @variable(model, [1:4bs, 1:4bs])
                            else
                                epos[i][j]= @variable(model, [1:bs, 1:bs])
                            end
                            for t = 1:bs, r = 1:bs
                                for s = 1:length(supp[k+1])
                                    a=ebasis[i][j][eblocks[i][j][t]]
                                    b=supp[k+1][s]
                                    c=ebasis[i][j][eblocks[i][j][r]]
                                    if ipart
                                        @inbounds bi = qtermadd3(b,c,a,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b,c,n)
                                    end
                                    bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                    Locb = bfind(tsupp, bi)
                                    if ipart == true
                                        Locb_cano = bfind(ttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])

                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][t,r])
                                    end
                                end
                                for s = 1:length(supp[k+1])
                                    a=ebasis[i][j][eblocks[i][j][t]]
                                    b=supp[k+1][s]
                                    c=ebasis[i][j][eblocks[i][j][r]]
                                    a_star = star(a,n)
                                    b_star = star(b,n)
                                    if ipart
                                        @inbounds bi = qtermadd3(c,a_star,b,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b_star,c,n)
                                    end
                                    bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                    Locb = bfind(tsupp, bi)
                                    if ipart == true
                                        Locb_cano = bfind(ttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])

                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])

                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][r,t])
                                    end
                                end
                            end
                        end
                    else
                        selected = ebasis[i][j][eblocks[i][j]]
                        bs = length(selected)
                        # Four independent real components for EVERY quaternion
                        # entry, including a singleton multiplier basis.
                        epos[i][j] = @variable(model, [1:4, 1:bs, 1:bs])
                        for t in 1:bs, r in 1:bs, s in eachindex(supp[k+1])
                            word = qtermadd3(supp[k+1][s], selected[r], selected[t], n)
                            h = ntuple(component -> epos[i][j][component,t,r], 4)
                            coefficient = quaternion_right_product(h, coe[k+1][s])
                            add_quaternion_equality_pair!(rcons, icons, jcons, kcons,
                                tsupp, ttsupp, word, coefficient, n, nb)
                        end
                    end
                end
                if addicons
                    addebasis = ebasis[i][end]
                    bs = length(addebasis)
                    if bs == 1
                        for rr = 1:3
                            for ss = 1:length(addsupp)
                                addepos[i][(rr-1)*length(addsupp)+ss] = @variable(model)
                                for tt= 1:4
                                    a = addebasis[1]
                                    b = addebasis[1]
                                    c = addsupp[ss][tt][1]
                                    d = addsupp[ss][tt][2]
                                    a_star = star(a,n)
                                    @inbounds bi = qtermadd4(d,b,a_star,c,n)
                                    bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                    Locb = bfind(tsupp, bi)
                                    if ipart
                                        Locb_cano = bfind(ttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                    end
                                end
                                for tt= 1:4
                                    a = addebasis[1]
                                    b = addebasis[1]
                                    c = addsupp[ss][tt][1]
                                    d = addsupp[ss][tt][2]
                                    a_star = star(a,n)
                                    c_star = star(c,n)
                                    d_star = star(d,n)
                                    @inbounds bi = qtermadd4(c_star,b,a_star,d_star,n)
                                    bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                    Locb = bfind(tsupp, bi)
                                    if ipart
                                        Locb_cano = bfind(ttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                    end
                                end
                            end
                        end
                    else
                        for rr = 1:3
                            for ss = 1:length(addsupp)
                                addepos[i][(rr-1)*length(addsupp)+ss] = @variable(model, [1:bs, 1:bs])
                                for t = 1:bs, r = 1:bs
                                    for tt= 1:4
                                        a = addebasis[t]
                                        b = addebasis[r]
                                        c = addsupp[ss][tt][1]
                                        d = addsupp[ss][tt][2]
                                        a_star = star(a,n)
                                        @inbounds bi = qtermadd4(d,b,a_star,c,n)
                                        bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                        Locb = bfind(tsupp, bi)
                                        if ipart == true
                                            Locb_cano = bfind(ttsupp, bi_cano)
                                            @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][t,r])

                                            @inbounds add_to_expression!(icons[Locb], imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss][t,r])

                                            @inbounds add_to_expression!(jcons[Locb], imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss][t,r])

                                            @inbounds add_to_expression!(kcons[Locb], imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss][t,r])

                                        else
                                            @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][t,r])
                                        end
                                    end
                                end
                                for t = 1:bs, r = 1:bs
                                    for tt = 1:4
                                        a = addebasis[t]
                                        b = addebasis[r]
                                        c = addsupp[ss][tt][1]
                                        d = addsupp[ss][tt][2]
                                        a_star = star(a,n)
                                        c_star = star(c,n)
                                        d_star = star(d,n)
                                        @inbounds bi = qtermadd4(c_star,b,a_star,d_star,n)
                                        bi, bi_cano = reduced_terms(bi, n, nb, ipart)
                                        Locb = bfind(tsupp, bi)
                                        if ipart == true
                                            Locb_cano = bfind(ttsupp, bi_cano)
                                            @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][r,t])

                                            @inbounds add_to_expression!(icons[Locb], -imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss][r,t])

                                            @inbounds add_to_expression!(jcons[Locb], -imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss][r,t])

                                            @inbounds add_to_expression!(kcons[Locb], -imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss][r,t])

                                        else
                                            @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][r,t])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        rbc = zeros(ipart ? lttsupp : ltsupp)
        if ipart == true
            ibc = zeros(ltsupp)
            jbc = zeros(ltsupp)
            kbc = zeros(ltsupp)
        end
        ncons = ltsupp
        if QUIET == false
            println("There are $ncons affine constraints.")
        end
        for i = 1:length(supp[1])
            mono = supp[1][i]
            mono_cano = canonical_qmono(mono,n)
            if ipart != true
                mono = canonical_qmono(mono,n)
            end

            if nb > 0
                mono = qreduce_unitnorm(mono,n, nb=nb)
                mono_cano = qreduce_unitnorm(mono_cano,n, nb=nb)
            end

            Locb = bfind(tsupp, mono)
            if ipart == true
                Locb_cano = bfind(ttsupp, mono_cano)
            end

            if Locb === nothing || (ipart == true && Locb_cano === nothing)
                throw(ArgumentError("The monomial basis is not enough for objective monomial $mono."))
            else
                if ipart == true
                    rbc[Locb_cano] += real(coe[1][i])
                    ibc[Locb] += imag_part(coe[1][i])[1]
                    jbc[Locb] += imag_part(coe[1][i])[2]
                    kbc[Locb] += imag_part(coe[1][i])[3]
                else
                    rbc[Locb] += real(coe[1][i])
                end
            end
        end
        @variable(model, lower)
        rcons[1] += lower
        inds = conjugate_representatives(tsupp, n; cyclic=!ipart)
        if ipart
            inds1 = conjugate_representatives(ttsupp, n; cyclic=true)
        end
        if ipart == true || qncnumeq > 0
            @constraint(model,rcon[k=1:length(inds1)],rcons[inds1[k]] == rbc[inds1[k]])
            @constraint(model, icon[k=1:length(inds)], icons[inds[k]] == ibc[inds[k]])
            @constraint(model, jcon[k=1:length(inds)], jcons[inds[k]] == jbc[inds[k]])
            @constraint(model, kcon[k=1:length(inds)], kcons[inds[k]] == kbc[inds[k]])
        else
            @constraint(model,rcon[k=1:length(inds)],rcons[inds[k]] == rbc[inds[k]])
        end

        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time = @elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        SDP_status = termination_status(model)
        objv = has_sdp_point(model) ? objective_value(model) : nothing
        if SDP_status != MOI.OPTIMAL
            println("termination status: $SDP_status")
            status = primal_status(model)
            println("solution status: $status")
        end
        QUIET || println("optimum = $objv")
    end
    return objv,ksupp,SDP_status
end

"""
    quaternion_to_real(pop, q)

Substitute real, i, j, and k coordinates for each quaternion and take the real
part. Return `(real_pop, x)` with coordinates grouped in four blocks of length
`length(q) ÷ 2`. The second half of `q` denotes conjugates of the first half.
"""
function quaternion_to_real(qpop,q)
    iseven(length(q)) && !isempty(q) || throw(ArgumentError("q must contain variables followed by their conjugates."))
    n = length(q) ÷ 2
    @polyvar x[1:4n]
    pop = Vector{Poly}(undef, length(qpop))
    I=quat(0,1,0,0)
    J=quat(0,0,1,0)
    K=quat(0,0,0,1)
    for (i,qp) in enumerate(qpop)
        if qp isa Number
            pop[i] = real(qp) + 0.0*x[1]
            continue
        end
        temp = qp(q[1:n]=>x[1:n]+I*x[n+1:2n]+J*x[2n+1:3n]+K*x[3n+1:4n], q[n+1:2n]=>x[1:n]-I*x[n+1:2n]-J*x[2n+1:3n]-K*x[3n+1:4n])
        pop[i] = real.(MultivariatePolynomials.coefficients(temp))'*MultivariatePolynomials.monomials(temp)
    end
    return pop,x
end

term_variables(term, n) = unique([v > n ? Int(v)-n : Int(v) for part in term for v in part])

function clique_decomp(n, m, rncnumeq, supp; alg="MF", qncnumeq=0, QUIET=false)
    if alg === false
        cliques = [collect(1:n)]
    else
        G = SimpleGraph(n)
        # Real commutation identities precede the final quaternion equalities.
        indices = [1:m+1-rncnumeq-qncnumeq; m+2-qncnumeq:m+1]
        for i in indices
            if i == 1
                for term in supp[i]
                    add_clique!(G, term_variables(term, n))
                end
            else
                # A multiplier belongs to a clique containing the full constraint.
                vars = unique([v for term in supp[i] for v in term_variables(term, n)])
                add_clique!(G, vars)
            end
        end
        cliques = alg == "NC" ? max_cliques(G)[1] : chordal_cliques!(G; method=alg, minimize=true)[1]
    end
    cliquesize = length.(cliques)
    QUIET || println("Variable clique sizes: $cliquesize")
    return cliques, length(cliques), cliquesize
end

function get_blocks(n, m, l, tsupp, supp::Vector{Vector{Vector{Vector{UInt16}}}}, basis, ebasis; nb=0, TS="block", merge=false, md=3)
    uk = m + 1
    blocks = Vector{Vector{Vector{Int}}}(undef, uk)
    blocksize = Vector{Vector{Int}}(undef, uk)
    cl = Vector{Int}(undef, uk)
    eblocks = Vector{Vector{Int}}(undef, l)
    if TS == false
        for k = 1:uk
            lb = length(basis[k])
            blocks[k],blocksize[k],cl[k] = [Vector(1:lb)],[lb],1
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else
        for k = 1:uk
            if k == 1
                G = get_graph(n, tsupp, basis[1], nb=nb)
            else
                G = get_graph(n, tsupp, supp[k-1], basis[k], nb=nb)
            end
            if TS == "block"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                if merge == true
                    blocks[k],cl[k],blocksize[k] = clique_merge!(blocks[k], d=md, QUIET=true)
                end
            end
        end
        for k = 1:l
            eblocks[k] = get_eblock(n, tsupp, supp[k+m], ebasis[k], nb=nb)
        end
    end
    return blocks,eblocks,cl,blocksize
end

function get_blocks(n, rlorder, I, J, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql, tsupp, basis, ebasis; TS="block", nb=0, merge=false, md=3)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    for i = 1:cql
        ksupp = TS == false ? nothing : tsupp[[
            let
                modified_3rd_l = [x > UInt16(n) ? x - UInt16(n) : x for x in tsupp[j][3]]
                issubset(union(tsupp[j][1], modified_3rd_l), cliques[i])
            end
            for j in eachindex(tsupp)
        ]]
        blocks[i],eblocks[i],cl[i],blocksize[i] = get_blocks(n, length(I[i]), length(J[i]), ksupp, supp[[I[i]; J[i]].+1], basis[i],
        ebasis[i], TS=TS, merge=merge, md=md, nb=nb)
    end
    return blocks,eblocks,cl,blocksize
end

function assign_constraint(n, m, numeq, supp, cliques, cql)
    I = [Int[] for _ in 1:cql]
    J = [Int[] for _ in 1:cql]
    ncc = Int[]
    for i in 1:m
        vars = unique([v > n ? Int(v)-n : Int(v) for term in supp[i+1] for part in term for v in part])
        ind = findall(k -> issubset(vars, cliques[k]), 1:cql)
        if isempty(ind)
            push!(ncc, i)
        elseif i <= m - numeq
            push!.(I[ind], i)
        else
            push!.(J[ind], i)
        end
    end
    return I, J, ncc
end

function get_graph(n, tsupp::Vector{Vector{Vector{UInt16}}}, basis; nb=0)
    lb = length(basis)
    G = SimpleGraph(lb)
    for i = 1:lb, j = 1:lb
        a=basis[i]
        b=basis[j]
        bi = qtermadd(a,b,n)
        if nb > 0
            bi = qreduce_unitnorm(bi, n, nb=nb)
        end
        if bfind(tsupp, bi) !== nothing
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_graph(n, tsupp::Vector{Vector{Vector{UInt16}}}, supp, basis; nb=0)
    lb = length(basis)
    G = SimpleGraph(lb)
    for i = 1:lb, j = 1:lb
        r = 1
        while r <= length(supp)
                a=basis[i]
                b=supp[r]
                c=basis[j]
                bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, n, nb=nb)
            end
            if bfind(tsupp, bi) !== nothing
                break
            else
                r += 1
            end
        end
        if r <= length(supp)
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_eblock(n, tsupp::Vector{Vector{Vector{UInt16}}}, hsupp::Vector{Vector{Vector{UInt16}}}, basis::Vector{Vector{Vector{UInt16}}}; nb=0)
    eblock = Int[]
    lb = length(basis)
    for i = 1:lb, j = 1:lb
        flag = 0
        for temp in hsupp
            a=basis[i]
            b=temp
            c=basis[j]
            bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, n, nb=nb)
            end
            if bfind(tsupp, bi) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            if i == j
                push!(eblock,i)
            else
                append!(eblock, [i,j])
            end
        end
    end
    return eblock
end

function get_gsupp(n, basis, ebasis, supp, cql, I, J, ncc, blocks, eblocks, cl, blocksize; nb=0)
    gsupp = Vector{Vector{UInt16}}[]
    for i = 1:cql
        for (j, w) in enumerate(I[i]), l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = 1:blocksize[i][j+1][l], item in supp[w+1]
            ind1 = blocks[i][j+1][l][t]
            ind2 = blocks[i][j+1][l][r]
            a=basis[i][j+1][ind1]
            b=item
            c=basis[i][j+1][ind2]
            @inbounds bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, n;nb=nb)
            end
            push!(gsupp, bi)
        end
        for (j, w) in enumerate(J[i]), k = 1: length(eblocks[i][j]),r = 1: length(eblocks[i][j]), item in supp[w+1]
            a=ebasis[i][j][eblocks[i][j][k]]
            b=item
            c=ebasis[i][j][eblocks[i][j][r]]
            @inbounds bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, n, nb=nb)
            end
            push!(gsupp, bi)
        end
    end
    for i in ncc, j = 1:length(supp[i+1])
        @inbounds bi = supp[i+1][j]
        if nb > 0
            bi = qreduce_unitnorm(bi, n, nb=nb)
        end
        push!(gsupp, bi)
    end
    return gsupp
end
