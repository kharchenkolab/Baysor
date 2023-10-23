using Statistics
using StatsBase

function estimate_hist(
        vec::Vector{<:Real}, weights=FrequencyWeights(ones(length(vec)));
        ext_cols::NamedTuple=NamedTuple(), rel_width::Float64=0.9, normalize::Union{Bool, Symbol}=false,
        center=false, bins=nothing, sigdigits::Int=4, type::Symbol=:rect, kwargs...
    )
    hf = (bins === nothing) ? fit(Histogram, vec, weights; kwargs...) : fit(Histogram, vec, weights, bins; kwargs...)
    diffs = rel_width * diff(hf.edges[1])[1]
    df = DataFrame(:s => hf.edges[1][1:end-1], :e => hf.edges[1][1:end-1] .+ diffs, :h => hf.weights)
    if center
        df.s .-= 0.5 * diffs
        df.e .-= 0.5 * diffs
    end
    for (k,v) in pairs(ext_cols)
        df[!,k] .= v
    end

    if isa(normalize, Bool)
        normalize = normalize ? :density : :none
    end

    if normalize == :density
        df[!, :h] = df.h ./ sum(df.h) ./ hf.edges[1].step.hi
    elseif normalize == :frac
        df[!, :h] = df.h ./ sum(df.h)
    elseif normalize != :none
        error("Unknown normalize")
    end

    if sigdigits > 0
        for c in [:s, :e, :h]
            df[!, c] = round.(df[!, c], sigdigits=sigdigits)
        end
    end

    if type == :rect
        df[!, :hs] .= 0.0
    elseif type == :bar
        select!(df, Not(:e))
    else
        error("Unknown type: $type. Only :rect and :bar are supported.")
    end

    return df
end

KWArgT = Union{Dict, NamedTuple, Nothing}

update_args(args::Union{Dict, NamedTuple}, update::Nothing) = args
update_args(args::Union{Dict, NamedTuple}, update::Union{Dict, NamedTuple}) =
    merge([Dict{Symbol, Any}(zip(keys(d), values(d))) for d in (args, update)]...)

function shuffle_labels(labels::Array{Int})
    new_labs = deepcopy(labels)
    mask = (new_labs .!= 0)
    new_labs[mask] = shuffle(1:maximum(labels))[new_labs[mask]]
    return new_labs
end

function shuffle_colors(colors::Vector)
    uniq_cols = unique(colors);
    col_ord = Dict(Pair.(uniq_cols, shuffle(1:length(uniq_cols))));
    return [uniq_cols[col_ord[tc]] for tc in colors]
end
