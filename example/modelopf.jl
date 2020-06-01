abstract type AbstractSparsePolyModel end

mutable struct SparsePolyModel <: AbstractSparsePolyModel
    n::UInt32
    m::UInt32
    numeq::UInt32
    supp::Vector{SparseMatrixCSC{UInt8,UInt32}}
    coe::Vector{Vector{Float64}}
    rlorder::Vector{Int}
end

function normalize(coe::Vector{Float64})
    mc = maximum(abs.(coe))
    return coe./mc
end

function move_zero(col,row,nz,coe)
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

function pop_opf(data::Dict{String, Any}; vmc="quartic")

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])

    nbus=length(ref[:bus])
    ng=length(ref[:gen])
    n=UInt32(2*nbus+2*ng)
#    m=UInt32(4*nbus+4*length(ref[:branch])+2*ng+length(ref[:ref_buses]))
    m=UInt32(3*nbus+4*length(ref[:branch])+2*ng+length(ref[:ref_buses]))
    numeq=UInt32(2*nbus+length(ref[:ref_buses]))
    rlorder=ones(Int, m+1)
    rlorder[1]=2
    supp=Vector{SparseMatrixCSC{UInt8,UInt32}}(undef, m+1)
    coe=Vector{Vector{Float64}}(undef, m+1)

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
    for (i,gen) in ref[:gen]
        coe[1][i+1]=gen["cost"][2]
        coe[1][i+ng+1]=gen["cost"][1]
    end
    col,row,nz,coe[1]=move_zero(col,row,nz,coe[1])
    supp[1]=SparseMatrixCSC(n,length(coe[1]),col,row,nz)

    bus=collect(keys(ref[:bus]))
    sort!(bus)
#    voltage magnitude constraints
    k=2
    if vmc=="quadratic"
        col=UInt32[1;1;2;3]
        nz=UInt8[2;2]
        for i=1:nbus
            row=UInt32[i;i+nbus]
            supp[k]=SparseMatrixCSC(n,3,col,row,nz)
            coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
            coe[k]=normalize(coe[k])
            k+=1
            supp[k]=SparseMatrixCSC(n,3,col,row,nz)
            coe[k]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
            coe[k]=normalize(coe[k])
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
            coe[k]=normalize(coe[k])
            rlorder[k]=0
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
        if maximum(bus)!=nbus
            vr_fr = bfind(bus,nbus,branch["f_bus"])
            vr_to = bfind(bus,nbus,branch["t_bus"])
        else
            vr_fr = branch["f_bus"]
            vr_to = branch["t_bus"]
        end
        vi_fr = vr_fr+nbus
        vi_to = vr_to+nbus

        # angle differences
        col=UInt32[1;3;5;7;9]
        nz=ones(UInt8,8)
        svr=sort([vr_fr;vr_to])
        svi=sort([vi_fr;vi_to])
        row=UInt32[svr;svi;vr_to;vi_fr;vr_fr;vi_to]
        coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
        coe[k]=normalize(coe[k])
        col,row,nz,coe[k]=move_zero(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1
        row=UInt32[vr_to;vi_fr;vr_fr;vi_to;svr;svi]
        coe[k]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
        coe[k]=normalize(coe[k])
        col,row,nz,coe[k]=move_zero(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        k+=1

        # thermal limits
        col=UInt32[1;1;2;4;6;8;11;13;16;18;21;24;26;28;29;31;33]
        row=UInt32[vr_fr;svr;vr_fr;vi_to;svr;svr;vi_fr;vr_fr;vi_fr;vr_fr;svi;vr_fr;vi_to;svr;vi_fr;vr_fr;svi;vr_to;vi_fr;vr_to;vi_fr;vi_fr;svi;svi]
        if vr_fr<vr_to&&vi_fr<vi_to
            nz=UInt8[4;3;1;3;1;2;2;2;1;1;2;2;2;1;1;2;2;1;1;2;1;2;1;2;2;1;3;4;3;1;2;2]
        elseif vr_fr<vr_to&&vi_fr>vi_to
            nz=UInt8[4;3;1;3;1;2;2;2;1;1;2;2;2;1;1;2;2;1;1;2;1;1;2;2;2;1;3;4;1;3;2;2]
        elseif vr_fr>vr_to&&vi_fr<vi_to
            nz=UInt8[4;1;3;3;1;2;2;1;2;1;2;2;2;1;1;2;2;1;1;2;1;2;1;2;2;1;3;4;3;1;2;2]
        else
            nz=UInt8[4;1;3;3;1;2;2;1;2;1;2;2;2;1;1;2;2;1;1;2;1;1;2;2;2;1;3;4;1;3;2;2]
        end
        coe[k]=[branch["rate_a"]^2*tm^4;ab1;acbd1;bcad1;cd1;-bcad1;2*ab1;acbd1;cd1;acbd1;bcad1;cd1;-bcad1;ab1;acbd1;cd1]
        coe[k]=normalize(coe[k])
        col,row,nz,coe[k]=move_zero(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        rlorder[k]=0
        k+=1
        col=UInt32[1;1;3;5;7;10;13;15;16;18;20;23;25;28;30;32;33]
        row=UInt32[svr;vr_fr;vi_to;svr;svr;vi_to;svr;vi_to;vr_fr;vi_to;vr_to;vr_to;vi_fr;vr_to;vi_fr;vr_to;svi;vr_to;vi_to;vr_to;svi;svi;svi;vi_to]
        if vr_fr<vr_to&&vi_fr<vi_to
            nz=UInt8[2;2;2;2;1;3;1;2;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;1;2;2;2;1;3;4]
        elseif vr_fr<vr_to&&vi_fr>vi_to
            nz=UInt8[2;2;2;2;1;3;1;2;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;2;1;2;2;3;1;4]
        elseif vr_fr>vr_to&&vi_fr<vi_to
            nz=UInt8[2;2;2;2;3;1;2;1;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;1;2;2;2;1;3;4]
        else
            nz=UInt8[2;2;2;2;3;1;2;1;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;2;1;2;2;3;1;4]
        end
        coe[k]=[branch["rate_a"]^2;cd2;cd2;acbd2;-bcad2;acbd2;-bcad2;ab2;bcad2;cd2;acbd2;2*ab2;bcad2;cd2;acbd2;ab2]
        coe[k]=normalize(coe[k])
        col,row,nz,coe[k]=move_zero(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        rlorder[k]=0
        k+=1
    end

    # power generation bound
    zero_pgen=UInt16[]
    zero_qgen=UInt16[]
    for (i, gen) in ref[:gen]
        if ref[:gen][i]["pmax"]>=1e-6
            col=UInt32[1;1;2;3]
            nz=UInt8[1;2]
            row=UInt32[i+2*nbus;i+2*nbus]
            coe[k]=[-ref[:gen][i]["pmin"]*ref[:gen][i]["pmax"];ref[:gen][i]["pmin"]+ref[:gen][i]["pmax"];-1]
            coe[k]=normalize(coe[k])
            col,row,nz,coe[k]=move_zero(col,row,nz,coe[k])
            supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
            k+=1
        else
            push!(zero_pgen, i)
            numeq+=1
        end
        if ref[:gen][i]["qmax"]>=1e-6
            col=UInt32[1;1;2;3]
            nz=UInt8[1;2]
            row=UInt32[i+2*nbus+ng;i+2*nbus+ng]
            coe[k]=[-ref[:gen][i]["qmin"]*ref[:gen][i]["qmax"];ref[:gen][i]["qmin"]+ref[:gen][i]["qmax"];-1]
            coe[k]=normalize(coe[k])
            col,row,nz,coe[k]=move_zero(col,row,nz,coe[k])
            supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
            k+=1
        else
            push!(zero_qgen, i)
            numeq+=1
        end
    end

    # active/reactive power
    for (i, bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        add_col=UInt32[2;4;6;8]
        add_nz=UInt8[1;1;1;1;1;1;1;1]
        col=UInt32[1;1;2;3]
        row=UInt32[i;i+nbus]
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
            if maximum(bus)!=nbus
                vr_fr = bfind(bus,nbus,branch["f_bus"])
                vr_to = bfind(bus,nbus,branch["t_bus"])
            else
                vr_fr = branch["f_bus"]
                vr_to = branch["t_bus"]
            end
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
            if vr_fr==i
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
                push!(bus_gen,gen_id)
            end
            lgen=length(bus_gen)
            append!(coe[k],-ones(Float64,lgen))
            append!(col,UInt32[i for i=1:lgen].+col[end])
            append!(row,bus_gen.+2*nbus)
            append!(nz,ones(UInt8,lgen))
            append!(coe[k+1],-ones(Float64,lgen))
            append!(qrow,bus_gen.+(2*nbus+ng))
        end
        qcol=copy(col)
        qnz=copy(nz)
        coe[k]=normalize(coe[k])
        col,row,nz,coe[k]=move_zero(col,row,nz,coe[k])
        supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
        coe[k+1]=normalize(coe[k+1])
        qcol,qrow,qnz,coe[k+1]=move_zero(qcol,qrow,qnz,coe[k+1])
        supp[k+1]=SparseMatrixCSC(n,length(coe[k+1]),qcol,qrow,qnz)
        k+=2
    end

    # reference voltage
    for key in keys(ref[:ref_buses])
        supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[key+nbus],UInt8[1])
        coe[k]=[1]
        k+=1
    end

    # zero power generation
    for i in zero_pgen
        supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+2*nbus],UInt8[1])
        coe[k]=[1]
        k+=1
    end
    for i in zero_qgen
        supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+2*nbus+ng],UInt8[1])
        coe[k]=[1]
        k+=1
    end
    return SparsePolyModel(n,m,numeq,supp,coe,rlorder)
end
