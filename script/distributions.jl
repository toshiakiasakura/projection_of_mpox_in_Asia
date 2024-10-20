using SpecialFunctions
using WeibullParetoDist
# Install via
# using Pkg
# Pkg.add(url="https://github.com/toshiakiasakura/WeibullParetoDist.jl")

include("util.jl")

NATSAL1Y_PARMS = [0.10, 0.77]
NATSAL4W_PARMS = [0.16, 0.87529]

rand_binom(n,p) = rand(Binomial(n,p))

"""Return scaled left truncated Weibull distributions.
...
# Arguments
- `scale_` : 365/infectious period should be given.
...
"""
function scaled_truncated_weibull(α, κ, scale_)
    κ_ = scale_^α*κ
    wb_trunc = truncated(WeibullPareto(α, κ_); lower= 1.0/scale_)
    return wb_trunc
end

function degree_and_probability_weibull(α::Float64, κ::Float64; k_max=10_000)
    wb = WeibullPareto(α, κ)
    wb_trunc = truncated(wb; lower=1.0)
    bins = LinRange( log(1), log(k_max), 51) .|> exp
    ks = [ mean(truncated(wb, bins[i], bins[i+1] )) for i in 1:(length(bins)-1)]
    Pk = [
        cdf(wb_trunc, bins[i+1]) - cdf(wb_trunc, bins[i])
        for i in 1:(length(bins)-1)
    ]
    return ks, Pk
end

"""Get an index pointing at the mean degree of the Weibull distribution
given one contact is taken randomly.
"""
function get_mean_index_for_simulation(α, κ; verbose=false, k_max=10000)
    wb = truncated(WeibullPareto(α, κ); lower=1.0)
    m = mean(wb)
    v = var(wb)
    m_degree = (v+m^2)/m

    bins = LinRange( log(1), log(k_max), 51) .|> exp
    ks, Pk = degree_and_probability_weibull(α, κ)
    cond = (bins .- m_degree) .< 0
    ind = sum(cond)
    if verbose == true
        println(m)
        println("Mean of p(x|T): ", m_degree)
        println("Fallen bins : ", bins[ind:ind+1])
        println("The correspoindg degree :",ks[ind])
    end
    return ind
end

"""
...
# Arguments
- `max_ind_k`: Number of compartment which population is larger than one
    given the whole population size.
    (Sometimes, high compartment does not have any).
- `k_max`: Maximum number of degree.
...
"""
function get_probs_given_a_contact(α, κ; max_ind_k=50, k_max=10_000)
    wb = truncated(WeibullPareto(α, κ); lower=1.0)
    m = mean(wb)
    v = var(wb)
    x = 1:10000
    k_max = 10_000

    bins = LinRange( log(1), log(k_max), 51) .|> exp
    target(x) = x*pdf(wb, x)/mean(wb)
    pmf = [quadgk(target, bins[i], bins[i+1])[1] for i in 1:(length(bins) -1)]
    if max_ind_k != 50
        pmf[max_ind_k+1:end] .= 0.0
    end
    return pmf
end