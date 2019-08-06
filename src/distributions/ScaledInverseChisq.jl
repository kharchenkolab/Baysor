using Distributions

struct ScaledInverseChisq <: Distributions.Distribution{Distributions.Univariate,Distributions.Continuous}
    ν::Float64
    τ::Float64
end

# Specify ScaledInverseChisq from mean and variance:
# μ = ντ² / (ν - 2)
# σ² = 2ν²τ⁴ / ((ν - 2)² (ν - 4))
# σ² = 2μ² / (ν - 4)
# ν = 2μ² / σ² + 4
# τ = √(μ (ν - 2) / ν)
function ScaledInverseChisq(; μ::Float64, σ²::Float64)
    ν = 2 * μ^2 / σ² + 4
    τ = sqrt(μ * (ν - 2) / ν)
    
    return ScaledInverseChisq(ν, τ)
end

rand(d::ScaledInverseChisq, args::Vararg{Int}) = d.τ^2 .* d.ν ./ rand(Chisq(d.ν), args...)
pdf(d::ScaledInverseChisq, x) = pdf(Chisq(d.ν), d.τ^2 * d.ν / x)