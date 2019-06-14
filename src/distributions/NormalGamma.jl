using Distributions
using StatsBase

"""
    NormalGamma(μ, κ, α, β)

Mean was estimated from κ observations with sample mean μ , and precision was estimated from 2α observations with sample mean μ and sum of squared deviations 2β.

# Arguments
- `μ::Float64`: mean of the distribution
- `κ::Float64`: number of observations, used to estimate mean
- `α::Float64`: 0.5 * number of observations used to estimate precision
- `β::Float64`: 0.5 * sum of squared deviations
"""
struct NormalGamma <: Distributions.Distribution{Distributions.Multivariate,Distributions.Continuous}
    μ::Float64;
    κ::Float64;
    α::Float64;
    β::Float64;

    t_scaled::TScaled;
    Γ_dist::Gamma;

    NormalGamma(μ::Real, κ::Real, α::Real, β::Real) =
        new(μ, κ, α, β, TScaled(μ, sqrt(β / (κ * α)), 2 * α), Gamma(α, 1 / β)) # Gamma here uses different parametrization
end

pdf(d::NormalGamma, μ::Float64, λ::Float64) = pdf(Normal(d.μ, 1 / sqrt(d.κ * λ)), μ) * pdf(d.Γ_dist, λ)
rand(d::NormalGamma) = (rand(d.t_scaled), rand(d.Γ_dist))
# var(d::NormalGamma) = d.β / d.α
