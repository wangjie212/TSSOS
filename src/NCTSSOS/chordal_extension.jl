function MCS(G)
    n=nv(G)
    w=zeros(UInt16,n)
    order=zeros(UInt16,n)
    for i=n:-1:1
        temp=w[[order[j]==0 for j=1:n]]
        v=findfirst(x->order[x]==0&&w[x]==maximum(temp), 1:n)
        w[neighbors(G, v)].+=1
        order[v]=i
    end
    return order
end

function MCSM!(G)
    n=nv(G)
    w=ones(UInt16,n)
    order=zeros(UInt16,n)
    for i=n:-1:1
        temp=w[[order[j]==0 for j=1:n]]
        maxw=maximum(temp)
        v=findfirst(x->order[x]==0&&w[x]==maxw, 1:n)
        order[v]=i
        unorder=findall(x->order[x]==0,1:n)
        unreach=copy(unorder)
        reach=[UInt16[] for j=1:maxw]
        S=UInt16[]
        for u in intersect(neighbors(G, v), unorder)
            push!(reach[w[u]], u)
            push!(S, u)
            deleteat!(unreach, bfind(unreach, length(unreach), u))
        end
        for k=1:maxw-1
            while reach[k]!=[]
                z=pop!(reach[k])
                for u in intersect(neighbors(G, z), unreach)
                    deleteat!(unreach, bfind(unreach, length(unreach), u))
                    if w[u]>k
                        push!(reach[w[u]], u)
                        push!(S, u)
                    else
                        push!(reach[k], u)
                    end
                end
            end
        end
        for u in S
            w[u]+=1
            add_edge!(G, v, u)
        end
    end
    return G,order
end

function Candidate(G,u,v,i,order,n)
    flag=true
    for z in neighbors(G, u)
        if order[z]>i&&has_edge(G, v, z)&&!has_edge(G, findfirst(x->order[x]==i, 1:n), z)
            flag=false
        end
    end
    return flag
end

function MinimalChordal!(G, order, F)
    n=nv(G)
    for i=n:-1:1
        ind=UInt16[]
        inc=UInt16[]
        for j=1:size(F[i],2)
            if Candidate(G,F[i][1,j],F[i][2,j],i,order,n)
                push!(ind, j)
                push!(inc, F[i][1,j], F[i][2,j])
            end
        end
        if length(ind)>0
            cand=F[i][:,ind]
            unique!(inc)
            sort!(inc)
            linc=length(inc)
            W=complete_graph(linc)
            for j=1:size(cand,2)
                rem_edge!(W, bfind(inc,linc,cand[1,j]), bfind(inc,linc,cand[2,j]))
            end
            keep,_ = MCSM!(W)
            for j=1:size(cand,2)
                if !has_edge(keep, bfind(inc,linc,cand[1,j]), bfind(inc,linc,cand[2,j]))
                    rem_edge!(G, cand[1,j], cand[2,j])
                end
            end
        end
    end
    return G,MCS(G)
end

function GreedyOrder!(G; method="MF", minimize=false)
    n=nv(G)
    H=copy(G)
    order=zeros(UInt16,n)
    if minimize==true
        F=Vector{Array{UInt16,2}}(undef, n)
    end
    for i=1:n
        if minimize==true
            F[i]=zeros(UInt16,2,1)
        end
        if method=="MF"
            num_fill=zeros(Int,n)
            for j=1:n
                neib=neighbors(H, j)
                lneib=length(neib)
                if lneib>0
                    sg=0
                    for k=1:lneib
                        for l=k+1:lneib
                            if has_edge(H, neib[k], neib[l])
                                sg+=1
                            end
                        end
                    end
                    num_fill[j]=Int(lneib*(lneib-1)/2-sg)
                elseif order[j]!=0
                    num_fill[j]=typemax(Int)
                end
            end
            v=argmin(num_fill)
        else
            deg=LightGraphs.degree(H)
            for j=1:n
                if order[j]!=0
                    deg[j]=typemax(Int)
                end
            end
            v=argmin(deg)
        end
        order[v]=i
        neib=copy(neighbors(H, v))
        for j=1:length(neib)
            rem_edge!(H, v, neib[j])
            for k=j+1:length(neib)
                if add_edge!(H, neib[j], neib[k])
                    add_edge!(G, neib[j], neib[k])
                    if minimize==true
                        F[i]=[F[i] [neib[j]; neib[k]]]
                    end
                end
            end
        end
        if minimize==true
            F[i]=F[i][:,2:end]
        end
    end
    if minimize==true
        G,order=MinimalChordal!(G, order, F)
    end
    return G,order
end

function chordal_cliques!(G; method="MF", minimize=false)
    G, order=GreedyOrder!(G, method=method, minimize=minimize)
    n=nv(G)
    candidate_cliques=Vector{Vector{UInt16}}(undef, n)
    for i=1:n
        inter=intersect(neighbors(G, i), findall(x->order[x]>order[i], 1:n))
        inter=[inter;i]
        candidate_cliques[i]=sort!(inter)
    end
    sort!(candidate_cliques, by = x -> length(x))
    reverse!(candidate_cliques)
    maximal_cliques = [first(candidate_cliques)]
    for clique in Iterators.drop(candidate_cliques, 1)
        if all(other_clique -> !(clique ⊆ other_clique), maximal_cliques)
            push!(maximal_cliques, clique)
        end
    end
    cliquesize=length.(maximal_cliques)
    cql=length(cliquesize)
    return maximal_cliques,cql,cliquesize
end

function add_clique!(G, nodes)
    for i in 1:length(nodes)-1
        for j in i+1:length(nodes)
            add_edge!(G, nodes[i], nodes[j])
        end
    end
end

function max_cliques(G)
    cliques=convert(Vector{Vector{UInt16}}, maximal_cliques(G))
    sort!.(cliques)
    cliquesize=length.(cliques)
    cql=length(cliquesize)
    return cliques,cql,cliquesize
end

# function fill!(G, order)
#     n=nv(G)
#     F=Vector{Array{UInt16,2}}(undef, n)
#     for i=1:n
#         v=findfirst(x->order[x]==i, 1:n)
#         neib=copy(intersect(neighbors(G, v), findall(x->order[x]>i, 1:n)))
#         F[i]=zeros(UInt16,2,1)
#         for j=1:length(neib)-1
#             for k=j+1:length(neib)
#                 if add_edge!(G, neib[j], neib[k])
#                     F[i]=[F[i] [neib[j]; neib[k]]]
#                 end
#             end
#         end
#         F[i]=F[i][:,2:end]
#     end
#     return G,F
# end
