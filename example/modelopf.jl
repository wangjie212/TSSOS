abstract type AbstractSparsePolyModel end

mutable struct SparsePolyModel <: AbstractSparsePolyModel
    n::UInt32
    m::UInt32
    numeq::UInt32
    nbus::Int
    ng::Int
    nb::Int
    supp::Vector{SparseMatrixCSC{UInt8,UInt32}}
    coe::Vector{Vector{Float64}}
    dg::Vector{Int}
end

function normalize(coe::Vector{Float64})
    mc = maximum(abs.(coe))
    return coe./mc
end

function move_zero!(col,row,nz,coe)
    i=1
    while i<=length(coe)
        if abs(coe[i])<=1e-8
            deleteat!(coe,i)
            deleteat!(row,col[i]:(col[i+1]-1))
            deleteat!(nz,col[i]:(col[i+1]-1))
            lrow=col[i+1]-col[i]
            deleteat!(col,i)
            col[i:end].-=lrow
        else
            i+=1
        end
    end
    col=convert(Vector{UInt32},col)
    return col,row,nz,coe
end

function fl_sum(vector)
    return mapreduce(x->x, +, vector, init = 0.0)
end

function pop_opf(data::Dict{String, Any}; vmc="quadratic",gen_model="two",normal=true)
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]
    nbus=length(ref[:bus])
    nb=length(ref[:branch])
    ng=length(ref[:gen])
    n=UInt32(2*nbus+2*ng)
    if vmc=="quadratic"&&gen_model=="two"
        m=UInt32(4*nbus+4*nb+2*ng+length(ref[:ref_buses]))
    elseif vmc=="quartic"&&gen_model=="two"
        m=UInt32(3*nbus+4*nb+2*ng+length(ref[:ref_buses]))
    elseif vmc=="quadratic"&&gen_model=="one"
        m=UInt32(4*nbus+4*nb+4*ng+length(ref[:ref_buses]))
    else
        m=UInt32(3*nbus+4*nb+4*ng+length(ref[:ref_buses]))
    end
    numeq=UInt32(2*nbus+length(ref[:ref_buses]))
    dg=2*ones(Int, m)
    supp=Vector{SparseMatrixCSC{UInt8,UInt32}}(undef, m+1)
    coe=Vector{Vector{Float64}}(undef, m+1)

    gens=collect(keys(ref[:gen]))
    sort!(gens)
    # objective function
    nc=2*ng+1
    col=UInt32[i for i=1:nc]
    col=[1;col]
    row=UInt32[i+2*nbus for i=1:ng]
    append!(row,row)
    nz=ones(UInt8,ng)
    append!(nz,2*nz)
    coe[1]=Vector{Float64}(undef, nc)
    coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
    for i=1:ng
        gen=ref[:gen][gens[i]]
        coe[1][i+1]=gen["cost"][2]
        coe[1][i+ng+1]=gen["cost"][1]
    end
    col,row,nz,coe[1]=move_zero!(col,row,nz,coe[1])
    supp[1]=SparseMatrixCSC(n,length(coe[1]),col,row,nz)

    bus=collect(keys(ref[:bus]))
    sort!(bus)
    # voltage magnitude constraints
    k=2
    if vmc=="quadratic"
        col=UInt32[1;1;2;3]
        nz=UInt8[2;2]
        for i=1:nbus
            row=UInt32[i;i+nbus]
            supp[k]=SparseMatrixCSC(n,3,col,row,nz)
            coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            k+=1
            supp[k]=SparseMatrixCSC(n,3,col,row,nz)
            coe[k]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            k+=1
        end
    else
        col=UInt32[1;1;2;3;4;6;7]
        nz=UInt8[2;2;4;2;2;4]
        for i=1:nbus
            row=UInt32[i;i+nbus;i;i;i+nbus;i+nbus]
            supp[k]=SparseMatrixCSC(n,6,col,row,nz)
            lb=ref[:bus][bus[i]]["vmin"]^2
            ub=ref[:bus][bus[i]]["vmax"]^2
            coe[k]=[-lb*ub;lb+ub;lb+ub;-1;-2;-1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            dg[k-1]=4
            k+=1
        end
    end

    for (i, branch) in ref[:branch]
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]
        ab1=-(g+g_fr)^2-(b+b_fr)^2
        cd1=-(-g*tr+b*ti)^2-(-b*tr-g*ti)^2
        acbd1=-2*((g+g_fr)*(-g*tr+b*ti)+(b+b_fr)*(-b*tr-g*ti))
        bcad1=-2*((b+b_fr)*(-g*tr+b*ti)-(g+g_fr)*(-b*tr-g*ti))
        ab2=-(g+g_to)^2-(b+b_to)^2
        cd2=-(-g*tr-b*ti)^2/tm^4-(-b*tr+g*ti)^2/tm^4
        acbd2=-2*((g+g_to)*(-g*tr-b*ti)+(b+b_to)*(-b*tr+g*ti))/tm^2
        bcad2=-2*((b+b_to)*(-g*tr-b*ti)-(g+g_to)*(-b*tr+g*ti))/tm^2
        vr_fr = bfind(bus,nbus,branch["f_bus"])
        vr_to = bfind(bus,nbus,branch["t_bus"])
        vi_fr = vr_fr+nbus
        vi_to = vr_to+nbus
        svr=sort([vr_fr;vr_to])
        svi=sort([vi_fr;vi_to])

        # angle differences
        col=UInt32[1;3;5;7;9]
        nz=ones(UInt8,8)
        row=UInt32[svr;svi;vr_to;vi_fr;vr_fr;vi_to]
        coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
        col=UInt32[1;3;5;7;9]
        nz=ones(UInt8,8)
        row=UInt32[vr_to;vi_fr;vr_fr;vi_to;svr;svi]
        coe[k]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1

        # thermal limits
        col=UInt32[1;1;2;4;6;8;11;13;16;18;21;24;26;28;29;31;33]
        row=UInt32[vr_fr;svr;vr_fr;vi_to;svr;svr;vi_fr;vr_fr;vi_fr;vr_fr;svi;vr_fr;vi_to;svr;vi_fr;vr_fr;svi;vr_to;vi_fr;vr_to;vi_fr;vi_fr;svi;svi]
        if vr_fr<vr_to
            nz=UInt8[4;3;1;3;1;2;2;2;1;1;2;2;2;1;1;2;2;1;1;2;1;2;1;2;2;1;3;4;3;1;2;2]
        else
            nz=UInt8[4;1;3;3;1;2;2;1;2;1;2;2;2;1;1;2;2;1;1;2;1;1;2;2;2;1;3;4;1;3;2;2]
        end
        coe[k]=[branch["rate_a"]^2*tm^4;ab1;acbd1;bcad1;cd1;-bcad1;2*ab1;acbd1;cd1;acbd1;bcad1;cd1;-bcad1;ab1;acbd1;cd1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        dg[k-1]=4
        k+=1
        col=UInt32[1;1;3;5;7;10;13;15;16;18;20;23;25;28;30;32;33]
        row=UInt32[svr;vr_fr;vi_to;svr;svr;vi_to;svr;vi_to;vr_fr;vi_to;vr_to;vr_to;vi_fr;vr_to;vi_fr;vr_to;svi;vr_to;vi_to;vr_to;svi;svi;svi;vi_to]
        if vr_fr<vr_to
            nz=UInt8[2;2;2;2;1;3;1;2;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;1;2;2;2;1;3;4]
        else
            nz=UInt8[2;2;2;2;3;1;2;1;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;2;1;2;2;3;1;4]
        end
        coe[k]=[branch["rate_a"]^2;cd2;cd2;acbd2;-bcad2;acbd2;-bcad2;ab2;bcad2;cd2;acbd2;2*ab2;bcad2;cd2;acbd2;ab2]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        dg[k-1]=4
        k+=1
   end

    # power generation bound
    zero_pgen=UInt16[]
    for i=1:ng
        gen=ref[:gen][gens[i]]
        if gen["pmax"]>=1e-6
            if gen_model=="two"
                col=UInt32[1;1;2;3]
                nz=UInt8[1;2]
                row=UInt32[i+2*nbus;i+2*nbus]
                coe[k]=[-gen["pmin"]*gen["pmax"];gen["pmin"]+gen["pmax"];-1]
                if normal==true
                    coe[k]=normalize(coe[k])
                end
                col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
                supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
                k+=1
            else
                col=UInt32[1;1;2]
                nz=UInt8[1]
                row=UInt32[i+2*nbus]
                coe[k]=[gen["pmax"];-1]
                if normal==true
                    coe[k]=normalize(coe[k])
                end
                supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
                k+=1
                col=UInt32[1;1;2]
                nz=UInt8[1]
                row=UInt32[i+2*nbus]
                coe[k]=[-gen["pmin"];1]
                if normal==true
                    coe[k]=normalize(coe[k])
                end
                col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
                supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
                k+=1
            end
        else
            push!(zero_pgen, i)
            numeq+=1
            if gen_model=="one"
                m-=1
            end
        end
        if gen_model=="two"
            col=UInt32[1;1;2;3]
            nz=UInt8[1;2]
            row=UInt32[i+2*nbus+ng;i+2*nbus+ng]
            coe[k]=[-gen["qmin"]*gen["qmax"];gen["qmin"]+gen["qmax"];-1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
            supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
            k+=1
        else
            col=UInt32[1;1;2]
            nz=UInt8[1]
            row=UInt32[i+2*nbus+ng]
            coe[k]=[gen["qmax"];-1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
            k+=1
            col=UInt32[1;1;2]
            nz=UInt8[1]
            row=UInt32[i+2*nbus]
            coe[k]=[-gen["qmin"];1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
            k+=1
        end
    end

    # active/reactive power
    for r=1:nbus
        i=bus[r]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        add_col=UInt32[2;4;6;8]
        add_nz=UInt8[1;1;1;1;1;1;1;1]
        col=UInt32[1;1;2;3]
        row=UInt32[r;r+nbus]
        nz=UInt8[2;2]
        coe[k]=zeros(Float64,4*length(ref[:bus_arcs][i])+3)
        coe[k+1]=zeros(Float64,4*length(ref[:bus_arcs][i])+3)
        coe[k][1]=fl_sum(load["pd"] for load in bus_loads)
        coe[k+1][1]=fl_sum(load["qd"] for load in bus_loads)
        sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)
        coe[k][2]=sgs
        coe[k][3]=sgs
        coe[k+1][2]=-sbs
        coe[k+1][3]=-sbs
        j=1
        for flow in ref[:bus_arcs][i]
            append!(col,col[end].+add_col)
            append!(nz,add_nz)
            branch=ref[:branch][flow[1]]
            vr_fr = bfind(bus,nbus,branch["f_bus"])
            vr_to = bfind(bus,nbus,branch["t_bus"])
            vi_fr = vr_fr+nbus
            vi_to = vr_to+nbus
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            temp1=(g+g_fr)/tm^2
            temp2=g+g_to
            temp3=(-g*tr+b*ti)/tm^2
            temp4=(-b*tr-g*ti)/tm^2
            temp5=(-g*tr-b*ti)/tm^2
            temp6=(-b*tr+g*ti)/tm^2
            temp7=-(b+b_fr)/tm^2
            temp8=-(b+b_to)
            if vr_fr==r
                coe[k][2]+=temp1
                coe[k][3]+=temp1
                coe[k][4+4*(j-1)]=temp3
                coe[k][5+4*(j-1)]=temp3
                coe[k][6+4*(j-1)]=temp4
                coe[k][7+4*(j-1)]=-temp4
                coe[k+1][2]+=temp7
                coe[k+1][3]+=temp7
                coe[k+1][4+4*(j-1)]=-temp4
                coe[k+1][5+4*(j-1)]=-temp4
                coe[k+1][6+4*(j-1)]=temp3
                coe[k+1][7+4*(j-1)]=-temp3
            else
                coe[k][2]+=temp2
                coe[k][3]+=temp2
                coe[k][4+4*(j-1)]=temp5
                coe[k][5+4*(j-1)]=temp5
                coe[k][6+4*(j-1)]=-temp6
                coe[k][7+4*(j-1)]=temp6
                coe[k+1][2]+=temp8
                coe[k+1][3]+=temp8
                coe[k+1][4+4*(j-1)]=-temp6
                coe[k+1][5+4*(j-1)]=-temp6
                coe[k+1][6+4*(j-1)]=-temp5
                coe[k+1][7+4*(j-1)]=temp5
            end
            if vr_fr<vr_to
                append!(row,[vr_fr;vr_to;vi_fr;vi_to;vr_to;vi_fr;vr_fr;vi_to])
            else
                append!(row,[vr_to;vr_fr;vi_to;vi_fr;vr_to;vi_fr;vr_fr;vi_to])
            end
            j+=1
        end
        qrow=copy(row)
        if !isempty(ref[:bus_gens][i])
            bus_gen=UInt32[]
            for gen_id in ref[:bus_gens][i]
                push!(bus_gen,bfind(gens,ng,gen_id))
            end
            lgen=length(bus_gen)
            append!(coe[k],-ones(Float64,lgen))
            append!(col,UInt32[l for l=1:lgen].+col[end])
            append!(row,bus_gen.+2*nbus)
            append!(nz,ones(UInt8,lgen))
            append!(coe[k+1],-ones(Float64,lgen))
            append!(qrow,bus_gen.+(2*nbus+ng))
        end
        qcol=copy(col)
        qnz=copy(nz)
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        if normal==true
            coe[k+1]=normalize(coe[k+1])
        end
        qcol,qrow,qnz,coe[k+1]=move_zero!(qcol,qrow,qnz,coe[k+1])
        supp[k+1]=SparseMatrixCSC(n,length(coe[k+1]),qcol,qrow,qnz)
        k+=2
    end

    # reference voltage
    for key in keys(ref[:ref_buses])
        i=bfind(bus,nbus,key)
        supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+nbus],UInt8[2])
        coe[k]=[1]
        k+=1
    end

    # zero power generation
    for i in zero_pgen
        supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+2*nbus],UInt8[1])
        coe[k]=[1]
        dg[k-1]=1
        k+=1
    end
    return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg)
end

function pop_opf_two(data::Dict{String, Any}; normal=true)
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]
    nbus=length(ref[:bus])
    ng=length(ref[:gen])
    nb=length(ref[:branch])
    n=UInt32(2*nbus+2*ng+4*nb)
    m=UInt32(4*nbus+8*length(ref[:branch])+2*ng+length(ref[:ref_buses]))
    numeq=UInt32(2*nbus+4*length(ref[:branch])+length(ref[:ref_buses]))
    dg=2*ones(Int, m)
    supp=Vector{SparseMatrixCSC{UInt8,UInt32}}(undef, m+1)
    coe=Vector{Vector{Float64}}(undef, m+1)

    gens=collect(keys(ref[:gen]))
    sort!(gens)
    # objective function
    nc=2*ng+1
    col=UInt32[i for i=1:nc]
    col=[1;col]
    row=UInt32[i+2*nbus for i=1:ng]
    append!(row,row)
    nz=ones(UInt8,ng)
    append!(nz,2*nz)
    coe[1]=Vector{Float64}(undef, nc)
    coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
    for i=1:ng
        gen=ref[:gen][gens[i]]
        coe[1][i+1]=gen["cost"][2]
        coe[1][i+ng+1]=gen["cost"][1]
    end
    col,row,nz,coe[1]=move_zero!(col,row,nz,coe[1])
    supp[1]=SparseMatrixCSC(n,length(coe[1]),col,row,nz)

    bus=collect(keys(ref[:bus]))
    sort!(bus)
    # voltage magnitude constraints
    k=2
    col=UInt32[1;1;2;3]
    nz=UInt8[2;2]
    for i=1:nbus
        row=UInt32[i;i+nbus]
        supp[k]=SparseMatrixCSC(n,3,col,row,nz)
        coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        k+=1
        supp[k]=SparseMatrixCSC(n,3,col,row,nz)
        coe[k]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        k+=1
    end

    branchs=collect(keys(ref[:branch]))
    sort!(branchs)
    for i=1:nb
        branch=ref[:branch][branchs[i]]
        vr_fr = bfind(bus,nbus,branch["f_bus"])
        vr_to = bfind(bus,nbus,branch["t_bus"])
        vi_fr = vr_fr+nbus
        vi_to = vr_to+nbus
        svr=sort([vr_fr;vr_to])
        svi=sort([vi_fr;vi_to])

        # angle differences
        col=UInt32[1;3;5;7;9]
        nz=ones(UInt8,8)
        row=UInt32[svr;svi;vr_to;vi_fr;vr_fr;vi_to]
        coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
        col=UInt32[1;3;5;7;9]
        nz=ones(UInt8,8)
        row=UInt32[vr_to;vi_fr;vr_fr;vi_to;svr;svi]
        coe[k]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1

        # thermal limits
        col=UInt32[1;1;2;3]
        row=UInt32[2*nbus+2*ng+i;2*nbus+2*ng+nb+i]
        nz=[2;2]
        coe[k]=[branch["rate_a"]^2;-1;-1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        supp[k]=SparseMatrixCSC(n,3,col,row,nz)
        coe[k+1]=coe[k]
        row=UInt32[2*nbus+2*ng+2*nb+i;2*nbus+2*ng+3*nb+i]
        supp[k+1]=SparseMatrixCSC(n,3,col,row,nz)
        k+=2
   end

    # power generation bound
    zero_pgen=UInt16[]
    for i=1:ng
        gen=ref[:gen][gens[i]]
        if gen["pmax"]>=1e-6
            col=UInt32[1;1;2;3]
            nz=UInt8[1;2]
            row=UInt32[i+2*nbus;i+2*nbus]
            coe[k]=[-gen["pmin"]*gen["pmax"];gen["pmin"]+gen["pmax"];-1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
            supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
            k+=1
        else
            push!(zero_pgen, i)
            numeq+=1
        end
        col=UInt32[1;1;2;3]
        nz=UInt8[1;2]
        row=UInt32[i+2*nbus+ng;i+2*nbus+ng]
        coe[k]=[-gen["qmin"]*gen["qmax"];gen["qmin"]+gen["qmax"];-1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
    end

    for i=1:nb
        branch=ref[:branch][branchs[i]]
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]
        vr_fr = bfind(bus,nbus,branch["f_bus"])
        vr_to = bfind(bus,nbus,branch["t_bus"])
        vi_fr = vr_fr+nbus
        vi_to = vr_to+nbus
        svr=sort([vr_fr;vr_to])
        svi=sort([vi_fr;vi_to])

        # Line Flow
        col=UInt32[1;2;3;5;7;9;11;12]
        nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
        row=UInt32[vr_fr;vi_fr;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+i]
        coe[k]=[g+g_fr;g+g_fr;-g*tr+b*ti;-g*tr+b*ti;-b*tr-g*ti;b*tr+g*ti;-tm^2]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
        col=UInt32[1;2;3;5;7;9;11;12]
        nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
        row=UInt32[vr_fr;vi_fr;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+nb+i]
        coe[k]=[-(b+b_fr);-(b+b_fr);b*tr+g*ti;b*tr+g*ti;-g*tr+b*ti;g*tr-b*ti;-tm^2]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
        col=UInt32[1;2;3;5;7;9;11;12]
        nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
        row=UInt32[vr_to;vi_to;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+2*nb+i]
        coe[k]=[(g+g_to)*tm^2;(g+g_to)*tm^2;-g*tr-b*ti;-g*tr-b*ti;b*tr-g*ti;-b*tr+g*ti;-tm^2]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
        col=UInt32[1;2;3;5;7;9;11;12]
        nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
        row=UInt32[vr_to;vi_to;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+3*nb+i]
        coe[k]=[-(b+b_to)*tm^2;-(b+b_to)*tm^2;b*tr-g*ti;b*tr-g*ti;g*tr+b*ti;-g*tr-b*ti;-tm^2]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
   end

    # active/reactive power
    for r=1:nbus
        i=bus[r]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        lg=length(ref[:bus_arcs][i])+length(ref[:bus_gens][i])+2
        col=UInt32[i for i=1:lg+1]
        col=UInt32[1;col]
        row=zeros(UInt32,lg)
        row[1]=r
        row[2]=r+nbus
        qrow=copy(row)
        nz=ones(UInt8,lg)
        nz[1]=2
        nz[2]=2
        coe[k]=ones(Float64,lg+1)
        coe[k+1]=ones(Float64,lg+1)
        coe[k][1]=fl_sum(load["pd"] for load in bus_loads)
        coe[k+1][1]=fl_sum(load["qd"] for load in bus_loads)
        sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)
        coe[k][2]=sgs
        coe[k][3]=sgs
        coe[k+1][2]=-sbs
        coe[k+1][3]=-sbs
        j=3
        for flow in ref[:bus_arcs][i]
            branch=ref[:branch][flow[1]]
            s=bfind(branchs,nb,flow[1])
            vr_fr = bfind(bus,nbus,branch["f_bus"])
            if vr_fr==r
                row[j]=2*nbus+2*ng+s
                qrow[j]=2*nbus+2*ng+nb+s
            else
                row[j]=2*nbus+2*ng+2*nb+s
                qrow[j]=2*nbus+2*ng+3*nb+s
            end
            j+=1
        end
        for gen_id in ref[:bus_gens][i]
            gen=bfind(gens,ng,gen_id)
            coe[k][j+1]=-1
            coe[k+1][j+1]=-1
            row[j]=gen+2*nbus
            qrow[j]=gen+2*nbus+ng
            j+=1
        end
        qcol=copy(col)
        qnz=copy(nz)
        if normal==true
            coe[k]=normalize(coe[k])
        end
        col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        if normal==true
            coe[k+1]=normalize(coe[k+1])
        end
        qcol,qrow,qnz,coe[k+1]=move_zero!(qcol,qrow,qnz,coe[k+1])
        supp[k+1]=SparseMatrixCSC(n,length(coe[k+1]),qcol,qrow,qnz)
        k+=2
    end

    # reference voltage
    for key in keys(ref[:ref_buses])
        i=bfind(bus,nbus,key)
        supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+nbus],UInt8[2])
        coe[k]=[1]
        k+=1
    end

    # zero power generation
    for i in zero_pgen
        supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+2*nbus],UInt8[1])
        coe[k]=[1]
        dg[k-1]=1
        k+=1
    end
    return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg)
end

function clique_opf_two(n,m,nbus,supp;alg="MF",minimize=true)
    G=SimpleGraph(n)
    for i=1:m+1
        for j = 1:supp[i].n
            add_clique!(G,supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)])
        end
    end
    for i=1:nbus
        add_edge!(G,i,i+nbus)
    end
    if alg=="NC"
        cliques,cql,cliquesize=max_cliques(G)
    else
        cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("The clique sizes of varibles:\n$uc\n$sizes")
    return cliques,cql,cliquesize
end

function clique_opf_four(m,nbus,nb,supp;vmc="quadratic",alg="MF",minimize=true)
    G=SimpleGraph(n)
    if vmc=="quadratic"
        t=2*nbus+4*nb
    else
        t=nbus+4*nb
    end
    for i=2:t+1
        for j = 1:supp[i].n
            add_clique!(G,supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)])
        end
    end
    if alg=="NC"
        cliques,cql,cliquesize=max_cliques(G)
    else
        cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
    end
    for i=m-numeq+1:m-numeq+2*nbus
        cql+=1
        rind=sort(unique(supp[i+1].rowval))
        push!(cliques, rind)
        push!(cliquesize, length(rind))
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("The clique sizes of varibles:\n$uc\n$sizes")
    return cliques,cql,cliquesize
end

# function approx_sol_opf(moment,n,cliques,cql,cliquesize)
#     qsol=Float64[]
#     lcq=sum(cliquesize)
#     A=zeros(lcq,n)
#     q=1
#     for k=1:cql
#         cqs=cliquesize[k]
#         if cqs==1
#             push!(qsol, moment[k])
#         else
#             F=eigen(moment[k], cqs:cqs)
#             temp=sqrt(F.values[1])*F.vectors[:,1]
#             for l=1:k-1
#                 inter=intersect(cliques[l], cliques[k])
#                 if inter!=[]
#                     flag=0
#                     for r in inter
#                         if l==1
#                             ind1=bfind(cliques[l],cliquesize[l],r)
#                         else
#                             ind1=sum(cliquesize[s] for s=1:l-1)+bfind(cliques[l],cliquesize[l],r)
#                         end
#                         ind2=bfind(cliques[k],cqs,r)
#                         if (qsol[ind1]>=1e-3&&temp[ind2]<=-1e-3)||(qsol[ind1]<=-1e-3&&temp[ind2]>=1e-3)
#                             temp=-temp
#                             flag=1
#                             break
#                         elseif (qsol[ind1]>=1e-3&&temp[ind2]>=1e-3)||(qsol[ind1]<=-1e-3&&temp[ind2]<=-1e-3)
#                             flag=1
#                             break
#                         end
#                     end
#                     if flag==1
#                         break
#                     end
#                 end
#             end
#             append!(qsol, temp)
#         end
#         for j=1:cqs
#             A[q,cliques[k][j]]=1
#             q+=1
#         end
#     end
#     return (A'*A)\(A'*qsol)
# end
