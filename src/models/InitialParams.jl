struct InitialParams
    centers::Array{Float64, 2}
    covs::Array{Array{Float64,2},1};
    assignment::Array{Int64,1};
    n_comps::Int;

    function InitialParams(centers::Array{T, 2} where T <: Real, covs::Union{Array{Array{Float64,2},1}, Real}, assignment::Array{Int64,1})
        if isa(covs, Real)
            covs = [[covs 0.; 0. covs] for i in 1:size(centers, 1)]
        else
            @assert size(centers, 1) == size(covs, 1)
        end

        return new(centers, adjust_cov_matrix.(covs), assignment, size(centers, 1))
    end
end