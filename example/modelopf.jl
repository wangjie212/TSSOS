using PowerModels
using TSSOS

abstract type AbstractSparsePolyModel end

mutable struct SparsePolyModel <: AbstractSparsePolyModel
    n
    m
    numeq
    nbus
    ng
    nb
    pop
    startpoint
end

function normalize(p::poly{Float64})
    return poly(p.supp, p.coe/maximum(abs.(p.coe)))
end

function normalize(p::cpoly{ComplexF64})
    return cpoly(p.supp, p.coe/maximum(abs.(p.coe)))
end

function move_zero(p::poly{Float64})
    ind = [abs(item) >= 1e-8 for item in p.coe]
    return poly(p.supp[ind], p.coe[ind])
end

function move_zero(p::cpoly{ComplexF64})
    ind = [abs(item) >= 1e-8 for item in p.coe]
    return cpoly(p.supp[ind], p.coe[ind])
end

function fl_sum(vector)
    return mapreduce(x->x, +, vector, init = 0.0)
end

# real POP formulization (voltage-power)
function pop_opf_real(data::Dict{String, Any}; normal=true, AngleCons=false, LineLimit=false)
    silence()
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][0]
    nbus = length(ref[:bus])
    nb = length(ref[:branch])
    ng = length(ref[:gen])
    n = 2*nbus + 2*ng
    m = 4*nbus + 2*ng + length(keys(ref[:ref_buses]))
    if AngleCons == true
        m += 2*nb
    end
    if LineLimit == true || LineLimit=="relax"
        m += 2*nb
    end
    numeq = 2*nbus + length(keys(ref[:ref_buses]))
    startpoint = zeros(n)
    pop = Vector{poly{Float64}}(undef, m+1)

    gens = collect(keys(ref[:gen]))
    sort!(gens)
    # objective function
    nc = 2*ng + 1
    coe = Vector{Float64}(undef, nc)
    supp = Vector{Vector{UInt16}}(undef, nc)
    coe[1] = sum(gen["cost"][3] for (i,gen) in ref[:gen])
    supp[1] = UInt16[]
    for i = 1:ng
        gen = ref[:gen][gens[i]]
        coe[2*(i-1)+2:2*(i-1)+3] = [gen["cost"][2]; gen["cost"][1]]
        supp[2*(i-1)+2:2*(i-1)+3] = [[2*nbus+i], [2*nbus+i;2*nbus+i]]
    end
    pop[1] = move_zero(poly(supp, coe))
    k = 2

    bus = collect(keys(ref[:bus]))
    sort!(bus)
    # voltage magnitude constraints
    for i = 1:nbus
        pop[k] = poly(Vector{UInt16}[[], [i;i], [i+nbus;i+nbus]], [-ref[:bus][bus[i]]["vmin"]^2;1;1])
        pop[k+1] = poly(Vector{UInt16}[[], [i;i], [i+nbus;i+nbus]], [ref[:bus][bus[i]]["vmax"]^2;-1;-1])
        startpoint[i] = sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
        startpoint[i+nbus] = sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
        k += 2
    end

    if AngleCons == true || LineLimit == true || LineLimit=="relax"
        for (i, branch) in ref[:branch]
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            vr = bfind(bus, branch["f_bus"])
            vt = bfind(bus, branch["t_bus"])
            srt = sort([vr;vt])
            dsrt = sort([vr;vr;vt;vt])
            ab1 = (g+g_fr)^2 + (b+b_fr)^2
            cd1 = (-g*tr+b*ti)^2 + (b*tr+g*ti)^2
            acbd1 = (g+g_fr)*(-g*tr+b*ti) - (b+b_fr)*(b*tr+g*ti)
            bcad1 = -(b+b_fr)*(-g*tr+b*ti) - (g+g_fr)*(b*tr+g*ti)
            ab2 = (g+g_to)^2*tm^4 + (b+b_to)^2*tm^4
            cd2 = (g*tr+b*ti)^2 + (-b*tr+g*ti)^2
            acbd2 = -(g+g_to)*tm^2*(g*tr+b*ti) + (b+b_to)*tm^2*(-b*tr+g*ti)
            bcad2 = (b+b_to)*tm^2*(g*tr+b*ti) + (g+g_to)*tm^2*(-b*tr+g*ti)

            # angle differences
            if AngleCons == true
                pop[k] = poly(Vector{UInt16}[srt, srt.+nbus, [vt;vr+nbus], [vr;vt+nbus]], [tan(branch["angmax"]);tan(branch["angmax"]);-1;1])
                pop[k+1] = poly(Vector{UInt16}[[vt;vr+nbus], [vr;vt+nbus], srt, srt.+nbus], [1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])])
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            end

            # thermal limits
            if LineLimit == true
                coe = [branch["rate_a"]^2*tm^4;-ab1;-2*ab1;-ab1;-cd1;-cd1;-cd1;-cd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*bcad1;2*bcad1;-2*bcad1;2*bcad1]
                supp = Vector{UInt16}[[], [vr;vr;vr;vr], [vr;vr;vr+nbus;vr+nbus], [vr+nbus;vr+nbus;vr+nbus;vr+nbus], [vt;vt;vr+nbus;vr+nbus], [vr;vr;vt+nbus;vt+nbus],
                dsrt, dsrt.+nbus, sort([vr;vr;vr;vt]), [vr;vr;sort([vr+nbus;vt+nbus])], [srt;vr+nbus;vr+nbus], sort([vr;vr;vr;vt]).+nbus, [sort([vr;vr;vt]);vr+nbus],
                [vr;vr;vr;vt+nbus], [vt;vr+nbus;vr+nbus;vr+nbus], [vr;sort([vr+nbus;vr+nbus;vt+nbus])]]
                pop[k] = move_zero(poly(supp, coe))
                coe = [branch["rate_a"]^2*tm^4;-ab2;-2*ab2;-ab2;-cd2;-cd2;-cd2;-cd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*bcad2;2*bcad2;-2*bcad2;2*bcad2]
                supp = Vector{UInt16}[[], [vt;vt;vt;vt], [vt;vt;vt+nbus;vt+nbus], [vt+nbus;vt+nbus;vt+nbus;vt+nbus], [vr;vr;vt+nbus;vt+nbus], [vt;vt;vr+nbus;vr+nbus],
                dsrt, dsrt.+nbus, sort([vt;vt;vt;vr]), [vt;vt;sort([vt+nbus;vr+nbus])], [srt;vt+nbus;vt+nbus], sort([vt;vt;vt;vr]).+nbus, [sort([vt;vt;vr]);vt+nbus],
                [vt;vt;vt;vr+nbus], [vr;vt+nbus;vt+nbus;vt+nbus], [vt;sort([vt+nbus;vt+nbus;vr+nbus])]]
                pop[k+1] = move_zero(poly(supp, coe))
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            elseif LineLimit == "relax"
                mvr = ref[:bus][bus[vr]]["vmin"]^2
                mvt = ref[:bus][bus[vt]]["vmin"]^2
                coe = [branch["rate_a"]^2*tm^4/mvr;-ab1;-ab1;-cd1;-cd1;-2*acbd1;2*bcad1;-2*bcad1;-2*acbd1]
                supp = Vector{UInt16}[[], [vr;vr], [vr+nbus;vr+nbus], [vt;vt], [vt+nbus;vt+nbus], srt, [vr;vt+nbus], [vt;vr+nbus], srt.+nbus]
                pop[k] = move_zero(poly(supp, coe))
                coe = [branch["rate_a"]^2*tm^4/mvt;-ab2;-ab2;-cd2;-cd2;-2*acbd2;2*bcad2;-2*bcad2;-2*acbd2]
                supp = Vector{UInt16}[[], [vt;vt], [vt+nbus;vt+nbus], [vr;vr], [vr+nbus;vr+nbus], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
                pop[k+1] = move_zero(poly(supp, coe))
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            end
        end
    end

    # power generation bound
    zero_pgen = UInt16[]
    for i = 1:ng
        gen = ref[:gen][gens[i]]
        if gen["pmax"] >= 1e-6
            coe = [-gen["pmin"]*gen["pmax"];gen["pmin"]+gen["pmax"];-1]
            supp = Vector{UInt16}[[], [i+2*nbus], [i+2*nbus;i+2*nbus]]
            pop[k] = move_zero(poly(supp, coe))
            if normal == true
                pop[k] = normalize(pop[k])
            end
            startpoint[i+2*nbus] = (gen["pmin"]+gen["pmax"])/2
            k += 1
        else
            push!(zero_pgen, i)
            startpoint[i+2*nbus] = 0
            numeq += 1
        end
        coe = [-gen["qmin"]*gen["qmax"];gen["qmin"]+gen["qmax"];-1]
        supp = Vector{UInt16}[[], [i+2*nbus+ng], [i+2*nbus+ng;i+2*nbus+ng]]
        pop[k] = move_zero(poly(supp, coe))
        if normal == true
            pop[k] = normalize(pop[k])
        end
        startpoint[i+2*nbus+ng] = 0
        k += 1
    end

    # active/reactive power
    for (r, i) in enumerate(bus)
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        coe1 = zeros(4*length(ref[:bus_arcs][i])+3)
        supp1 = Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
        supp1[1:3] = [[], [r;r], [r+nbus;r+nbus]]
        coe1[1] = fl_sum(load["pd"] for load in bus_loads)
        coe2 = zeros(4*length(ref[:bus_arcs][i])+3)
        supp2 = Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
        supp2[1:3] = [[], [r;r], [r+nbus;r+nbus]]
        coe2[1] = fl_sum(load["qd"] for load in bus_loads)
        sgs = fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs = fl_sum(shunt["bs"] for shunt in bus_shunts)
        coe1[2:3] = [sgs;sgs]
        coe2[2:3] = [-sbs;-sbs]
        j = 1
        for flow in ref[:bus_arcs][i]
            branch = ref[:branch][flow[1]]
            vr = bfind(bus, branch["f_bus"])
            vt = bfind(bus, branch["t_bus"])
            srt = sort([vr;vt])
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            a1 = (g+g_fr)/tm^2
            b1 = -(b+b_fr)/tm^2
            c1 = (-g*tr+b*ti)/tm^2
            d1 = (b*tr+g*ti)/tm^2
            a2 = g + g_to
            b2 = -(b+b_to)
            c2 = -(g*tr+b*ti)/tm^2
            d2 = -(-b*tr+g*ti)/tm^2

            supp1[j+3:j+6] = [srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
            supp2[j+3:j+6] = [srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
            if vr == r
                coe1[2:3] .+= a1
                coe1[j+3:j+6] = [c1;-d1;d1;c1]
                coe2[2:3] .+= b1
                coe2[j+3:j+6] = [d1;c1;-c1;d1]
            else
                coe1[2:3] .+= a2
                coe1[j+3:j+6] = [c2;d2;-d2;c2]
                coe2[2:3] .+= b2
                coe2[j+3:j+6] = [d2;-c2;c2;d2]
            end
            j += 4
        end
        if !isempty(ref[:bus_gens][i])
            for gen_id in ref[:bus_gens][i]
                gen = bfind(gens, gen_id)
                push!(supp1, [gen+2*nbus])
                push!(coe1, -1)
                push!(supp2, [gen+2*nbus+ng])
                push!(coe2, -1)
            end
        end
        pop[k] = move_zero(poly(supp1, coe1))
        pop[k+1] = move_zero(poly(supp2, coe2))
         if normal == true
            pop[k] = normalize(pop[k])
            pop[k+1] = normalize(pop[k+1])
        end
        k += 2
    end

    # reference voltage
    for key in keys(ref[:ref_buses])
        i = bfind(bus, key)
        pop[k] = poly(Vector{UInt16}[[i+nbus;i+nbus]], Float64[1])
        k += 1
    end

    # zero power generation
    for i in zero_pgen
        pop[k] = poly(Vector{UInt16}[[i+2*nbus]], Float64[1])
        k += 1
    end

    return SparsePolyModel(n, m, numeq, nbus, ng, nb, pop, startpoint)
end

# complex POP formulization (voltage only)
function pop_opf_com_vol(data::Dict{String, Any}; normal=true, AngleCons=false, LineLimit=false)
    silence()
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][0]
    nbus = length(ref[:bus])
    nb = length(ref[:branch])
    ng = length(ref[:gen])
    n = nbus
    m = 4*nbus + 2*ng + length(keys(ref[:ref_buses]))
    if AngleCons == true
        m += 2*nb
    end
    if LineLimit == true || LineLimit == "relax"
        m += 2*nb
    end
    numeq = 2*nbus - 2*ng + length(keys(ref[:ref_buses]))
    pop = Vector{cpoly{ComplexF64}}(undef, m+1)

    gens = collect(keys(ref[:gen]))
    sort!(gens)
    # objective function
    pop[1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], [])], ComplexF64[sum(gen["cost"][3] for (i,gen) in ref[:gen])])
    k = 2

    bus = collect(keys(ref[:bus]))
    sort!(bus)
    # voltage magnitude constraints
    for i = 1:nbus
        pop[k] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([i], [i])], ComplexF64[-ref[:bus][bus[i]]["vmin"]^2;1])
        pop[k+1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([i], [i])], ComplexF64[ref[:bus][bus[i]]["vmax"]^2;-1])
        k += 2
    end

    if AngleCons == true || LineLimit == true || LineLimit == "relax"
        for (i, branch) in ref[:branch]
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            vr = bfind(bus, branch["f_bus"])
            vt = bfind(bus, branch["t_bus"])
            srt = sort([vr;vt])
            a1 = g + g_fr
            b1 = -(b+b_fr)
            c1 = -g*tr + b*ti
            d1 = b*tr + g*ti
            a2 = (g+g_to)*tm^2
            b2 = -(b+b_to)*tm^2
            c2 = -(g*tr+b*ti)
            d2 = -(-b*tr+g*ti)

            # angle differences
            if AngleCons == true
                pop[k] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([vr], [vt]), tuple([vt], [vr])], [tan(branch["angmax"])+im;tan(branch["angmax"])-im])
                pop[k+1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([vr], [vt]), tuple([vt], [vr])], [-tan(branch["angmin"])-im;-tan(branch["angmin"])+im])
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            end

            # thermal limits
            if LineLimit == true
                coe = [branch["rate_a"]^2*tm^4;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vr;vr], [vr;vr]), tuple([vr;vr], srt), tuple(srt, [vr;vr]), tuple(srt, srt)]
                pop[k] = move_zero(cpoly(supp, coe))
                coe = [branch["rate_a"]^2*tm^4;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vt;vt], [vt;vt]), tuple(srt, [vt;vt]), tuple([vt;vt], srt), tuple(srt, srt)]
                pop[k+1] = move_zero(cpoly(supp, coe))
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            elseif LineLimit == "relax"
                mvr = ref[:bus][bus[vr]]["vmin"]^2
                mvt = ref[:bus][bus[vt]]["vmin"]^2
                coe = [branch["rate_a"]^2*tm^4/mvr;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vr], [vr]), tuple([vr], [vt]), tuple([vt], [vr]), tuple([vt], [vt])]
                pop[k] = move_zero(cpoly(supp, coe))
                coe = [branch["rate_a"]^2*tm^4/mvt;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vt], [vt]), tuple([vr], [vt]), tuple([vt], [vr]), tuple([vr], [vr])]
                pop[k+1] = move_zero(cpoly(supp, coe))
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            end
        end
    end

    # active/reactive power
    k1 = k
    k = k1 + 4*ng
    for (r, i) in enumerate(bus)
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        ploads = fl_sum(load["pd"] for load in bus_loads)
        qloads = fl_sum(load["qd"] for load in bus_loads)
        sgs = fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs = fl_sum(shunt["bs"] for shunt in bus_shunts)
        if length(ref[:bus_gens][i]) > 1
            @error "There are more than one generator at one bus!"
        elseif isempty(ref[:bus_gens][i])
            coe1 = zeros(ComplexF64, 2*length(ref[:bus_arcs][i])+2)
            supp1 = Vector{Tuple{Vector{UInt16}, Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
            supp1[1:2] = [tuple([], []), tuple([r], [r])]
            coe1[1:2] = [ploads;sgs]
            coe2 = zeros(ComplexF64, 2*length(ref[:bus_arcs][i])+2)
            supp2 = Vector{Tuple{Vector{UInt16}, Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
            supp2[1:2] = [tuple([], []), tuple([r], [r])]
            coe2[1:2] = [qloads;-sbs]
        else
            gen_id = ref[:bus_gens][i][1]
            gen = ref[:gen][gen_id]
            gc2 = gen["cost"][2]
            pop[k1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([r], [r])], ComplexF64[gen["pmax"]-ploads, -sgs])
            pop[k1+1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([r], [r])], ComplexF64[-gen["pmin"]+ploads, sgs])
            pop[k1+2] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([r], [r])], ComplexF64[gen["qmax"]-qloads, sbs])
            pop[k1+3] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([r], [r])], ComplexF64[-gen["qmin"]+qloads, -sbs])
            push!(pop[1].supp, tuple([r], [r]))
            ngen = length(pop[1].supp)
            pop[1].coe[1] += gc2*ploads
            push!(pop[1].coe, gc2*sgs)
        end
        j = 1
        for flow in ref[:bus_arcs][i]
            branch = ref[:branch][flow[1]]
            vr = bfind(bus, branch["f_bus"])
            vt = bfind(bus, branch["t_bus"])
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            a1 = (g+g_fr)/tm^2
            b1 = -(b+b_fr)/tm^2
            c1 = (-g*tr+b*ti)/tm^2
            d1 = (b*tr+g*ti)/tm^2
            a2 = g + g_to
            b2 = -(b+b_to)
            c2 = -(g*tr+b*ti)/tm^2
            d2 = -(-b*tr+g*ti)/tm^2

            if isempty(ref[:bus_gens][i])
                supp1[j+2:j+3] = [tuple([vr], [vt]), tuple([vt], [vr])]
                supp2[j+2:j+3] = [tuple([vr], [vt]), tuple([vt], [vr])]
            else
                push!(pop[1].supp, tuple([vr], [vt]), tuple([vt], [vr]))
                push!(pop[k1].supp, tuple([vr], [vt]), tuple([vt], [vr]))
                push!(pop[k1+1].supp, tuple([vr], [vt]), tuple([vt], [vr]))
                push!(pop[k1+2].supp, tuple([vr], [vt]), tuple([vt], [vr]))
                push!(pop[k1+3].supp, tuple([vr], [vt]), tuple([vt], [vr]))
            end
            if vr == r
                if isempty(ref[:bus_gens][i])
                    coe1[2] += a1
                    coe1[j+2:j+3] = [(c1+d1*im)/2, (c1-d1*im)/2]
                    coe2[2] += b1
                    coe2[j+2:j+3] = [(-c1*im+d1)/2, (c1*im+d1)/2]
                else
                    pop[1].coe[ngen] += gc2*a1
                    push!(pop[1].coe, gc2*(c1+d1*im)/2, gc2*(c1-d1*im)/2)
                    pop[k1].coe[2] -= a1
                    push!(pop[k1].coe, -(c1+d1*im)/2, -(c1-d1*im)/2)
                    pop[k1+1].coe[2] += a1
                    push!(pop[k1+1].coe, (c1+d1*im)/2, (c1-d1*im)/2)
                    pop[k1+2].coe[2] -= b1
                    push!(pop[k1+2].coe, -(-c1*im+d1)/2, -(c1*im+d1)/2)
                    pop[k1+3].coe[2] += b1
                    push!(pop[k1+3].coe, (-c1*im+d1)/2, (c1*im+d1)/2)
                end
            else
                if isempty(ref[:bus_gens][i])
                    coe1[2] += a2
                    coe1[j+2:j+3] = [(c2-d2*im)/2, (c2+d2*im)/2]
                    coe2[2] += b2
                    coe2[j+2:j+3] = [(c2*im+d2)/2, (-c2*im+d2)/2]
                else
                    pop[1].coe[ngen] += gc2*a2
                    push!(pop[1].coe, gc2*(c2-d2*im)/2, gc2*(c2+d2*im)/2)
                    pop[k1].coe[2] -= a2
                    push!(pop[k1].coe, -(c2-d2*im)/2, -(c2+d2*im)/2)
                    pop[k1+1].coe[2] += a2
                    push!(pop[k1+1].coe, (c2-d2*im)/2, (c2+d2*im)/2)
                    pop[k1+2].coe[2] -= b2
                    push!(pop[k1+2].coe, -(c2*im+d2)/2, -(-c2*im+d2)/2)
                    pop[k1+3].coe[2] += b2
                    push!(pop[k1+3].coe, (c2*im+d2)/2, (-c2*im+d2)/2)
                end
            end
            j += 2
        end
        if isempty(ref[:bus_gens][i])
            pop[k] = move_zero(cpoly(supp1, coe1))
            pop[k+1] = move_zero(cpoly(supp2, coe2))
            if normal == true
                pop[k] = normalize(pop[k])
                pop[k+1] = normalize(pop[k+1])
            end
            k += 2
        else
            pop[k1] = move_zero(pop[k1])
            pop[k1+1] = move_zero(pop[k1+1])
            pop[k1+2] = move_zero(pop[k1+2])
            pop[k1+3] = move_zero(pop[k1+3])
            if normal == true
                pop[k1] = normalize(pop[k1])
                pop[k1+1] = normalize(pop[k1+1])
                pop[k1+2] = normalize(pop[k1+2])
                pop[k1+3] = normalize(pop[k1+3])
            end
            k1 += 4
            if gen["cost"][1] > 0
                lsupp = length(pop[1].supp)
                gc1 = gen["cost"][1]
                pop[1].coe[1] += gc1*ploads^2
                for l = ngen:lsupp
                    push!(pop[1].supp, tuple([pop[1].supp[l][1];pop[1].supp[l][1]], [pop[1].supp[l][2];pop[1].supp[l][2]]))
                    push!(pop[1].coe, gc1*pop[1].coe[l]^2/gc2^2)
                    for p = l+1:lsupp
                        push!(pop[1].supp, tuple(sort([pop[1].supp[l][1];pop[1].supp[p][1]]), sort([pop[1].supp[l][2];pop[1].supp[p][2]])))
                        push!(pop[1].coe, 2*gc1*pop[1].coe[l]*pop[1].coe[p]/gc2^2)
                    end
                    pop[1].coe[l] += 2*gc1*ploads*pop[1].coe[l]/gc2
                end
            end
        end
    end
    pop[1] = arrange(move_zero(pop[1]))

    # reference voltage
    for key in keys(ref[:ref_buses])
        i = bfind(bus, key)
        pop[k] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([i;i], []), tuple([i], [i]), tuple([], [i;i])], ComplexF64[1;-2;1])
        k += 1
    end
    return SparsePolyModel(n, m, numeq, nbus, ng, nb, pop, nothing)
end

# complex POP formulization (voltage-power)
function pop_opf_com(data::Dict{String, Any}; normal=true, AngleCons=false, LineLimit=false)
    silence()
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][0]
    nbus = length(ref[:bus])
    nb = length(ref[:branch])
    ng = length(ref[:gen])
    n = nbus + ng
    m = 4*nbus + 2*ng + length(keys(ref[:ref_buses]))
    if AngleCons == true
        m += 2*nb
    end
    if LineLimit == true || LineLimit=="relax"
        m += 2*nb
    end
    numeq = 2*nbus + length(keys(ref[:ref_buses]))
    pop = Vector{cpoly{ComplexF64}}(undef, m+1)

    gens = collect(keys(ref[:gen]))
    sort!(gens)
    # objective function
    nc = 5*ng + 1
    coe = Vector{ComplexF64}(undef, nc)
    supp = Vector{Tuple{Vector{UInt16}, Vector{UInt16}}}(undef, nc)
    coe[1] = sum(gen["cost"][3] for (i,gen) in ref[:gen])
    supp[1] = tuple([], [])
    for i = 1:ng
        gen = ref[:gen][gens[i]]
        coe[5*(i-1)+2:5*(i-1)+3] = [0.5*gen["cost"][2], 0.5*gen["cost"][2]]
        coe[5*(i-1)+4:5*(i-1)+6] = [0.25*gen["cost"][1], 0.5*gen["cost"][1], 0.25*gen["cost"][1]]
        supp[5*(i-1)+2:5*(i-1)+6] = [tuple([nbus+i], []), tuple([], [nbus+i]), tuple([nbus+i;nbus+i], []), tuple([nbus+i], [nbus+i]), tuple([], [nbus+i;nbus+i])]
    end
    pop[1] = move_zero(cpoly(supp, coe))
    k = 2

    bus = collect(keys(ref[:bus]))
    sort!(bus)
    # voltage magnitude constraints
    for i = 1:nbus
        pop[k] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([i], [i])], ComplexF64[-ref[:bus][bus[i]]["vmin"]^2;1])
        pop[k+1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([i], [i])], ComplexF64[ref[:bus][bus[i]]["vmax"]^2;-1])
        k += 2
    end

    if AngleCons == true || LineLimit == true || LineLimit == "relax"
        for (i, branch) in ref[:branch]
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            vr = bfind(bus, branch["f_bus"])
            vt = bfind(bus, branch["t_bus"])
            srt = sort([vr, vt])
            a1 = g + g_fr
            b1 = -(b+b_fr)
            c1 = - g*tr + b*ti
            d1 = b*tr + g*ti
            a2 = (g+g_to)*tm^2
            b2 = -(b+b_to)*tm^2
            c2 = -(g*tr+b*ti)
            d2 = -(-b*tr+g*ti)

            # angle differences
            if AngleCons == true
                pop[k] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([vr], [vt]), tuple([vt], [vr])], [tan(branch["angmax"])+im;tan(branch["angmax"])-im])
                pop[k+1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([vr], [vt]), tuple([vt], [vr])], [-tan(branch["angmin"])-im;-tan(branch["angmin"])+im])
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            end

            # thermal limits
            if LineLimit == true
                coe = [branch["rate_a"]^2*tm^4;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vr;vr], [vr;vr]), tuple([vr;vr], srt), tuple(srt, [vr;vr]), tuple(srt, srt)]
                pop[k] = move_zero(cpoly(supp, coe))
                coe = [branch["rate_a"]^2*tm^4;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vt;vt], [vt;vt]), tuple(srt, [vt;vt]), tuple([vt;vt], srt), tuple(srt, srt)]
                pop[k+1] = move_zero(cpoly(supp, coe))
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            elseif LineLimit == "relax"
                mvr = ref[:bus][bus[vr]]["vmin"]^2
                mvt = ref[:bus][bus[vt]]["vmin"]^2
                coe = [branch["rate_a"]^2*tm^4/mvr;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vr], [vr]), tuple([vr], [vt]), tuple([vt], [vr]), tuple([vt], [vt])]
                pop[k] = move_zero(cpoly(supp, coe))
                coe = [branch["rate_a"]^2*tm^4/mvt;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
                supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([vt], [vt]), tuple([vr], [vt]), tuple([vt], [vr]), tuple([vr], [vr])]
                pop[k+1] = move_zero(cpoly(supp, coe))
                if normal == true
                    pop[k] = normalize(pop[k])
                    pop[k+1] = normalize(pop[k+1])
                end
                k += 2
            end
        end
    end

    # power generation bound
    zero_pgen = UInt16[]
    for i = 1:ng
        gen = ref[:gen][gens[i]]
        if gen["pmax"] >= 1e-6
            coe = ComplexF64[-4*gen["pmin"]*gen["pmax"];2*gen["pmin"]+2*gen["pmax"];2*gen["pmin"]+2*gen["pmax"];-1;-2;-1]
            supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([i+nbus], []), tuple([], [i+nbus]), tuple([i+nbus;i+nbus], []), tuple([i+nbus], [i+nbus]), tuple([], [i+nbus;i+nbus])]
            pop[k] = move_zero(cpoly(supp, coe))
            if normal == true
                pop[k] = normalize(pop[k])
            end
            k += 1
        else
            push!(zero_pgen, i)
            numeq += 1
        end
        coe = [-4*gen["qmin"]*gen["qmax"];-2*gen["qmin"]*im-2*gen["qmax"]*im;2*gen["qmin"]*im+2*gen["qmax"]*im;1;-2;1]
        supp = Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([], []), tuple([i+nbus], []), tuple([], [i+nbus]), tuple([i+nbus;i+nbus], []), tuple([i+nbus], [i+nbus]), tuple([], [i+nbus;i+nbus])]
        pop[k] = move_zero(cpoly(supp, coe))
        if normal == true
            pop[k] = normalize(pop[k])
        end
        k += 1
    end

    # active/reactive power
    for (r, i) in enumerate(bus)
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        coe1 = zeros(ComplexF64, 2*length(ref[:bus_arcs][i])+2)
        supp1 = Vector{Tuple{Vector{UInt16}, Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
        supp1[1:2] = [tuple([], []), tuple([r], [r])]
        coe1[1] = fl_sum(load["pd"] for load in bus_loads)
        coe1[2] = fl_sum(shunt["gs"] for shunt in bus_shunts)
        coe2 = zeros(ComplexF64, 2*length(ref[:bus_arcs][i])+2)
        supp2 = Vector{Tuple{Vector{UInt16}, Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
        supp2[1:2] = [tuple([], []), tuple([r], [r])]
        coe2[1] = fl_sum(load["qd"] for load in bus_loads)
        coe2[2] = -fl_sum(shunt["bs"] for shunt in bus_shunts)
        j = 1
        for flow in ref[:bus_arcs][i]
            branch = ref[:branch][flow[1]]
            vr = bfind(bus, branch["f_bus"])
            vt = bfind(bus, branch["t_bus"])
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            a1 = (g+g_fr)/tm^2
            b1 = -(b+b_fr)/tm^2
            c1 = (-g*tr+b*ti)/tm^2
            d1 = (b*tr+g*ti)/tm^2
            a2 = g + g_to
            b2 = -(b+b_to)
            c2 = -(g*tr+b*ti)/tm^2
            d2 = -(-b*tr+g*ti)/tm^2

            supp1[j+2:j+3] = [tuple([vr], [vt]), tuple([vt], [vr])]
            supp2[j+2:j+3] = [tuple([vr], [vt]), tuple([vt], [vr])]
            if vr == r
                coe1[2] += a1
                coe1[j+2:j+3] = [(c1+d1*im)/2, (c1-d1*im)/2]
                coe2[2] += b1
                coe2[j+2:j+3] = [(-c1*im+d1)/2, (c1*im+d1)/2]
            else
                coe1[2] += a2
                coe1[j+2:j+3] = [(c2-d2*im)/2, (c2+d2*im)/2]
                coe2[2] += b2
                coe2[j+2:j+3] = [(c2*im+d2)/2, (-c2*im+d2)/2]
            end
            j += 2
        end
        if !isempty(ref[:bus_gens][i])
            for gen_id in ref[:bus_gens][i]
                gen = bfind(gens, gen_id)
                push!(supp1, tuple([gen+nbus], []), tuple([], [gen+nbus]))
                push!(coe1, -0.5, -0.5)
                push!(supp2, tuple([gen+nbus], []), tuple([], [gen+nbus]))
                push!(coe2, 0.5*im, -0.5*im)
            end
        end
        pop[k] = move_zero(cpoly(supp1, coe1))
        pop[k+1] = move_zero(cpoly(supp2, coe2))
        if normal == true
            pop[k] = normalize(pop[k])
            pop[k+1] = normalize(pop[k+1])
        end
        k += 2
    end

    # reference voltage
    for key in keys(ref[:ref_buses])
        i = bfind(bus, key)
        pop[k] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([i;i], []), tuple([i], [i]), tuple([], [i;i])], ComplexF64[1;-2;1])
        k += 1
    end

    # zero power generation
    for i in zero_pgen
        pop[k] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[tuple([i+nbus], []), tuple([], [i+nbus])], ComplexF64[1;1])
        k += 1
    end

    return SparsePolyModel(n, m, numeq, nbus, ng, nb, pop, nothing)
end
