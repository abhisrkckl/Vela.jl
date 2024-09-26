export KINPriorDistribution,
    SINIPriorDistribution, STIGMAPriorDistribution, SHAPMAXPriorDistribution

abstract type CustomPriorDistribution <: ContinuousUnivariateDistribution end

pdf(d::CustomPriorDistribution, x::Real) = insupport(d, x) ? pdf_expr(d, x) : 0.0
logpdf(d::CustomPriorDistribution, x::Real) = log(pdf(d, x))

function cdf(d::CustomPriorDistribution, x::Real)
    return if x <= minimum(d)
        0.0
    elseif x >= maximum(d)
        1.0
    else
        cdf_expr(d, x)
    end
end
logcdf(d::CustomPriorDistribution, x::Real) = log(cdf(d, x))

minimum(d::CustomPriorDistribution) = support(d).lb
maximum(d::CustomPriorDistribution) = support(d).ub


"""
    KINPriorDistribution

Distribution of KIN = ι when cos(ι) is uniformly distributed in [0,1].
"""
struct KINPriorDistribution <: CustomPriorDistribution end
support(::KINPriorDistribution) = RealInterval(0.0, π / 2)
pdf_expr(::KINPriorDistribution, ι) = sin(ι)
cdf_expr(::KINPriorDistribution, ι) = 1 - cos(ι)
quantile(::KINPriorDistribution, q) = acos(1 - q)


"""
    SINIPriorDistribution

Distribution of SINI = sin(ι) when cos(ι) is uniformly distributed in [0,1].
"""
struct SINIPriorDistribution <: CustomPriorDistribution end
support(::SINIPriorDistribution) = RealInterval(0.0, 1.0)
pdf_expr(::SINIPriorDistribution, s) = s / sqrt(1 - s * s)
cdf_expr(::SINIPriorDistribution, s) = 1 - sqrt(1 - s * s)
quantile(::SINIPriorDistribution, q) = sqrt((2 - q) * q)


"""
    STIGMAPriorDistribution

Distribution of STIGMA = sin(ι)/(1 + cos(ι)) when cos(ι) is uniformly distributed in [0,1].
"""
struct STIGMAPriorDistribution <: CustomPriorDistribution end
support(::STIGMAPriorDistribution) = RealInterval(0.0, 1.0)
pdf_expr(::STIGMAPriorDistribution, Ϛ) = 4 * Ϛ / (1 + Ϛ^2)^2
cdf_expr(::STIGMAPriorDistribution, Ϛ) = 2 * Ϛ^2 / (1 + Ϛ^2)
quantile(::STIGMAPriorDistribution, q) = q / (2 - q)


"""
    SHAPMAXPriorDistribution

Distribution of SHAPMAX = -ln(1 - sin(ι)) when cos(ι) is uniformly distributed in [0,1].
"""
struct SHAPMAXPriorDistribution <: CustomPriorDistribution end
support(::SHAPMAXPriorDistribution) = RealInterval(0.0, Inf)
pdf_expr(::SHAPMAXPriorDistribution, S) = (1 - exp(-S)) / (sqrt(2 * exp(S) - 1))
cdf_expr(::SHAPMAXPriorDistribution, S) = 1 - exp(-S) * sqrt(2 * exp(S) - 1)
quantile(::SHAPMAXPriorDistribution, q) = log((sqrt(2q - q^2) + 1) / (q - 1)^2)
