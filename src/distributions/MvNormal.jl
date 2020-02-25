using Distributions
using LinearAlgebra
using PDMats: PDMat, PDiagMat
using StatsBase

MvNormalF = MvNormal{Float64, T, Vector{Float64}} where T <: Union{PDMat{Float64, Matrix{Float64}}, PDiagMat{Float64, Vector{Float64}}}

function log_pdf(μx::Float64, μy::Float64, Σ::PDiagMat{Float64, Vector{Float64}}, x::Float64, y::Float64)::Float64
    std_x = sqrt(Σ.diag[1])
    std_y = sqrt(Σ.diag[2])

    ndx = (x - μx) / std_x
    ndy = (y - μy) / std_y

    log_div2 = log(2 * π * std_x * std_y)

    return -(ndx * ndx + ndy * ndy) / 2 - log_div2
end

@inbounds function log_pdf(μx::Float64, μy::Float64, Σ::Array{Float64,2}, x::Float64, y::Float64)::Float64
    std_x = sqrt(Σ[1, 1])
    std_y = sqrt(Σ[2, 2])
    ρ = max(min(1 - 1e-2, Σ[1, 2] / (std_x * std_y)), -1 + 1e-2)

    ndx = (x - μx) / std_x
    ndy = (y - μy) / std_y;

    div1 = 2 * (1 - ρ^2)
    log_div2 = log(2 * π * std_x * std_y * sqrt(1 - ρ^2))

    return -(ndx * ndx + ndy * ndy - 2 * ρ * ndx * ndy) / div1 - log_div2
end

log_pdf(μx::Float64, μy::Float64, Σ::PDMat{Float64,Array{Float64,2}}, x::Float64, y::Float64) = log_pdf(μx, μy, Σ.mat, x, y)
log_pdf(d::MvNormalF, x::Float64, y::Float64)::Float64 = log_pdf(d.μ[1]::Float64, d.μ[2]::Float64, d.Σ, x, y)

pdf(d::MvNormalF, x::Float64, y::Float64) = exp(log_pdf(d, x, y))
shape(d::MvNormalF) = Matrix(d.Σ)

function robust_cov(x::Array{Float64, 2}; prop::Float64=0.1)
    v1 = mad(x[:,1], normalize=true)^2
    v2 = mad(x[:,2], normalize=true)^2

    dev1 = x[:,1] .- trim_mean(x[:,1]; prop=prop)
    dev2 = x[:,2] .- trim_mean(x[:,2]; prop=prop)
    covar = mean(winsor(dev1 .* dev2; prop=prop))

    return [v1 covar; covar v2]
end

function maximize(dist::MvNormalF, x::Array{Float64, 2})::MvNormalF  # TODO: replace with maximize!
    if size(x, 2) == 0
        return dist
    end

    x = x'
    # μ = [trim_mean(x[:,1], prop=0.2), trim_mean(x[:,2], prop=0.2)]
    μ = vec(mean(x, dims=1))
    if size(x, 1) <= 2
        return dist
    end

    # Σ = adjust_cov_matrix!(robust_cov(x; prop=0.1))
    Σ = adjust_cov_matrix!(cov(x))

    if det(Σ) ≈ 0
        error("Singular covariance matrix")
    end

    return MvNormal(μ, Σ)
end

function adjust_cov_matrix!(Σ::Array{Float64, 2}; max_iter::Int=100,
                            cov_modifier::Float64=1e-4, tol=1e-10)::Array{Float64, 2}
    if isposdef(Σ) && det(Σ) > tol
        return Σ
    end

    for it in 1:max_iter
        if isposdef(Σ) && det(Σ) > tol
            break
        end

        @inbounds Σ[diagind(Σ)] .+= cov_modifier
        cov_modifier *= 1.5
    end

    return Σ
end