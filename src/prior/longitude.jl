export UniformLongitude

import Distributions: ContinuousUnivariateDistribution, @distr_support, pdf, logpdf, quantile

struct UniformLongitude <: ContinuousUnivariateDistribution end

pdf(::UniformLongitude, x::Real) = (-π/2 < x < π/2) ?  (cos(x) / 2) : 0.0
logpdf(ul::UniformLongitude, x::Real) = log(pdf(ul, x))

# function cdf(::UniformLongitude, x::Real)
#     if x < -π/2
#         return 0.0
#     elseif x > π/2
#         return 1.0
#     else
#         return (1 + sin(x)) / 2
#     end
# end

quantile(::UniformLongitude, q::Real) = asin(2*q - 1)
