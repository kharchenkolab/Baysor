using Distributions
using LinearAlgebra
using PDMats
using StatsBase

function log_pdf(μx::Float64, μy::Float64, Σ::PDMats.PDiagMat{Float64,Array{Float64,1}}, x::Float64, y::Float64)::Float64
    std_x = sqrt(Σ.diag[1])
    std_y = sqrt(Σ.diag[2])

    ndx = (x - μx) / std_x
    ndy = (y - μy) / std_y

    log_div2 = log(2 * π * std_x * std_y)

    return -(ndx * ndx + ndy * ndy) / 2 - log_div2
end

function log_pdf(μx::Float64, μy::Float64, Σ::Array{Float64,2}, x::Float64, y::Float64)::Float64
    std_x = sqrt(Σ[1, 1])
    std_y = sqrt(Σ[2, 2])
    ρ = max(min(1 - 1e-2, Σ[1, 2] / (std_x * std_y)), -1 + 1e-2)

    ndx = (x - μx) / std_x
    ndy = (y - μy) / std_y;

    div1 = 2 * (1 - ρ^2)
    log_div2 = log(2 * π * std_x * std_y * sqrt(1 - ρ^2))

    return -(ndx * ndx + ndy * ndy - 2 * ρ * ndx * ndy) / div1 - log_div2
end

log_pdf(μx::Float64, μy::Float64, Σ::PDMats.PDMat{Float64,Array{Float64,2}}, x::Float64, y::Float64) = log_pdf(μx, μy, Σ.mat, x, y)
log_pdf(d::MvNormal, x::Float64, y::Float64)::Float64 = log_pdf(d.μ[1], d.μ[2], d.Σ, x, y)

pdf(d::MvNormal, x::Float64, y::Float64) = exp(log_pdf(d, x, y))
shape(d::MvNormal) = Matrix(d.Σ)

function robust_cov(x::Array{Float64, 2}; prop::Float64=0.1)
    v1 = mad(x[:,1], normalize=true)^2
    v2 = mad(x[:,2], normalize=true)^2

    dev1 = x[:,1] .- trim_mean(x[:,1]; prop=prop)
    dev2 = x[:,2] .- trim_mean(x[:,2]; prop=prop)
    covar = mean(winsor(dev1 .* dev2; prop=prop))

    return [v1 covar; covar v2]
end

function maximize(dist::Distributions.MvNormal, x::Array{Float64, 2})::Distributions.MvNormal
    if size(x, 2) == 0
        return dist
    end

    x = copy(x')
    μ = [trim_mean(x[:,1], prop=0.2), trim_mean(x[:,2], prop=0.2)]
    if size(x, 1) <= 2
        return dist
    end

    Σ = adjust_cov_matrix(robust_cov(x; prop=0.1))

    if det(Σ) ≈ 0
        @show Σ
        throw(ErrorException("Singular covariance matrix"))
    end

    return MvNormal(μ, Σ)
end

function adjust_cov_matrix(Σ::Array{Float64, 2}; max_iter::Int=100,
                           cov_modifier::Array{Float64, 2}=Matrix(Diagonal(1e-4 .* ones(size(Σ, 2)))), tol=1e-10)::Array{Float64, 2}
    if isposdef(Σ) && det(Σ) > tol
        return Σ
    end

    res = copy(Σ)
    for it in 1:max_iter
        if isposdef(res) && det(res) > tol
            break
        end

        res += cov_modifier
        cov_modifier *= 1.5
    end

    return res
end