struct InitialParams{L}
    centers::Matrix{Float64}
    covs::Array{CovMat{L},1};
    assignment::Vector{Int64};

    InitialParams(centers::Matrix{<:Real}, covs::Float64, assignment::Vector{Int64}) =
        new{N}(deepcopy(centers), [CovMat{N}(diagm(0 => (ones(N) .* covs))) for i in 1:size(centers, 1)], deepcopy(assignment))

    InitialParams(centers::Matrix{<:Real}, covs::Array{<:CovMat{N}, 1}, assignment::Array{Int64,1}) where N =
        new{N}(deepcopy(centers), deepcopy(covs), deepcopy(assignment))
end