using Distributions

struct TScaled <: Distributions.Distribution{Distributions.Univariate,Distributions.Continuous}
    μ::Float64;
    σ::Float64;
    ν::Float64;

    t_dist::TDist;

    TScaled(μ::Real, σ::Real, ν::Real) = new(μ, σ, ν, TDist(ν))
end

pdf(d::TScaled, x::Float64) = pdf(d.t_dist, (x - d.μ) / d.σ)
rand(d::TScaled) = rand(d.t_dist) * d.σ + d.μ
