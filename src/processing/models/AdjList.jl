import Base: getindex, length, eachindex, keys

struct AdjList
    ids::Vector{Vector{Int}}
    weights::Vector{Vector{Float64}}
end

getindex(adj_list::AdjList, i::Int) = (adj_list.ids[i], adj_list.weights[i])
length(adj_list::AdjList) = length(adj_list.ids)
eachindex(adj_list::AdjList) = eachindex(adj_list.ids)
keys(adj_list::AdjList) = keys(adj_list.ids)