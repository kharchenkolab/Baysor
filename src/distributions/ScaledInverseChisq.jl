using Distributions

struct ScaledInverseChisq <: Distributions.Distribution{Distributions.Univariate,Distributions.Continuous}
    ν::Int
    τ::Float64
end

rand(d::ScaledInverseChisq, args::Vararg{Int}) = d.τ^2 .* d.ν ./ rand(Chisq(d.ν), args...)
pdf(d::ScaledInverseChisq, x) = pdf(Chisq(d.ν), d.τ^2 * d.ν / x)