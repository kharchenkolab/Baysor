using SparseArrays

struct InitialParams
    centers::Matrix{Float64}
    covs::Array{CovMat,1};
    assignment::Vector{Int64};

    InitialParams(centers::Array{T, 2} where T <: Real, covs::Float64, assignment::Array{Int64,1}) =
        new(deepcopy(centers), [CovMat([covs 0.; 0. covs]) for i in 1:size(centers, 1)], deepcopy(assignment))

    InitialParams(centers::Array{T, 2} where T <: Real, covs::Array{CovMat, 1}, assignment::Array{Int64,1}) =
        new(deepcopy(centers), deepcopy(covs), deepcopy(assignment))
end