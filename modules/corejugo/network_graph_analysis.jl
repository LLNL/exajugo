# function to build network graph

function build_network_graph(sys::SCACOPFdata,
                             contingency::Union{GenericContingency, Nothing}=nothing)
    
    # create simple graph with all buses as vertices
    g = SimpleGraph(size(sys.N, 1))
    
    # add all lines as edges
    for l = 1:size(sys.L, 1)
        if !isnothing(contingency) && l in contingency.lines_out
            continue
        end
        add_edge!(g, sys.L_Nidx[l, 1], sys.L_Nidx[l, 2])
    end
    
    # add all transformers as edges
    for t = 1:size(sys.T, 1)
        if !isnothing(contingency) && t in contingency.transformers_out
            continue
        end
        add_edge!(g, sys.T_Nidx[t, 1], sys.T_Nidx[t, 2])
    end
    
    # return graph object
    return g
    
end

# function to determine which elements of a vector are in a second 

velems_in_velems(x::Vector, y::Vector)::BitVector = .!isnothing.(indexin(x, y))

# function to get bus indices affected by contingencies

function contingency_bus_indices(sys::SCACOPFdata)::Vector{Int}
    out = Vector{Int}(undef, size(sys.K, 1))
    for i = 1:size(sys.K, 1)
        if sys.K[i, :ConType] == :Generator
            out[i] = sys.G_Nidx[sys.K[i, :IDout]]
        elseif sys.K[i, :ConType] == :Line
            out[i] = sys.L_Nidx[sys.K[i, :IDout], 1]
        else
            @assert sys.K[i, :ConType] == :Transformer
            out[i] = sys.T_Nidx[sys.K[i, :IDout], 1]
        end
    end
    return out
end

# function to get number of subsystems

function number_of_connected_subsystems(sys::SCACOPFdata)::Int
    g = build_network_graph(sys)
    cc = connected_components(g)
    return length(cc)
end

# function to separate system into connected subsystems

function split_on_connected_subsystems(sys::SCACOPFdata)::Vector{SCACOPFdata}
    
    # create network graph
    g = build_network_graph(sys)
    
    # get connected components
    cc = connected_components(g)
    
    # if there is only 1 connected component, we have nothing to do
    if length(cc) == 1
        println("System has 1 connected component. No splitting necessary.")
        return [sys]
    end
    
    # create a new SCACOPFdata for each connected component
    out = SCACOPFdata[]
    con_bus_idx = contingency_bus_indices(sys)
    for cci in cc
        Ncci = sys.N[cci, :]
        Lcci = sys.L[velems_in_velems(sys.L_Nidx[:, 1], cci), :]
        Tcci = sys.T[velems_in_velems(sys.T_Nidx[:, 1], cci), :]
        SShcci = sys.SSh[velems_in_velems(sys.SSh_Nidx, cci), :]
        Gcci = sys.G[velems_in_velems(sys.G_Nidx, cci), :]
        Kcci = sys.K[velems_in_velems(con_bus_idx, cci), :]
        sys_cci = SCACOPFdata(sys.MVAbase, Ncci, Lcci, Tcci, SShcci, Gcci, Kcci,
                              sys.P, sys.generators, DataFrame([sys.cont_labels], [:LABEL]))
        push!(out, sys_cci)
    end
    
    # print splitting summary
    println("Splitted system onto ", length(cc),
            " subsystems, with total bus counts from ", minimum(length.(cc)),
            " to ", maximum(length.(cc)), ".")
    
    # return vector with connected component systems
    return out
    
end

