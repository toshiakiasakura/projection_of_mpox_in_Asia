# -*- coding: utf-8 -*-
# ## Adjust assortativity based on reported Newman's assortativity coefficient

include("./utils.jl")
include("./distributions.jl")
include("./model.jl")
include("./model_intercountry.jl")

KS_4W, PK_4W = degree_and_probability_weibull(NATSAL4W_PARMS...)
scale_v = 1/4
println("Number of sexual activity group with less than 100 sexual partnership per 3 months: ", sum((KS_4W * scale_v) .< 100))

# Check the number of sexual partnerships per 3 months to match EricChow2016. 
inds = [(1, 11), (12, 15), (16, 33), (34, 50)]
for ind in inds
    ks = round.([KS_4W[ind[1]], KS_4W[ind[2]]] .* scale_v, digits=1)
    println("Bins: $(ks[1]) ~ $(ks[2]) ")
end

KS_4W, PK_4W = degree_and_probability_weibull(NATSAL4W_PARMS...)
normalize!(PK_4W, 1)
mat_e50 = normalize(KS_4W .* PK_4W, 1) .* reshape(normalize(KS_4W .* PK_4W,1), 1, 50)
@test round(sum(mat_e50) , digits=5) == 1
mat_e3 = convert_50comp_to_EricChow2016_bins(mat_e50, KS_4W)
Newman_assortativity_coef(mat_e3)

heatmap(mat_e50, right_margin=5Plots.mm, 
    xlabel="Sexual activity group of contactor", 
    ylabel="Sexual activity group of contactee",
    fmt=:png,
)

# Version 1: Only increase a monoganous partner. 
α = 1.8
mat_e50_α = increase_low_assortativity(mat_e50, α, 1.)
mat_e3_α = convert_50comp_to_EricChow2016_bins(mat_e50_α, KS_4W)
println("Newman_assortativity_coef: ", Newman_assortativity_coef(mat_e3_α))

rem_prop = 0.3015 * (α - 1)/α
println("Proportion of MSM population which should be removed from a simulation (%): ", rem_prop*100, " and reduced to be: ", 1- rem_prop)

heatmap(mat_e50_α, fmt=:png,
    xlabel="Sexual activity group of contactor", 
    ylabel="Sexual activity group of contactee",
)

# Version 2: Increase both 1-1 and 2-3 to 2-3. 
α1 = 2.00
α2_3 = 5.41
mat_e50_α = increase_low_assortativity(mat_e50, α1, α2_3)
mat_e3_α = convert_50comp_to_EricChow2016_bins(mat_e50_α, KS_4W)
Newman_assortativity_coef(mat_e3_α)

# +
# Proportion of rewired total partnerhsips among all partnerships
mat_e50_α = increase_low_assortativity(mat_e50, α1, α2_3; norm=false)
calculate_proportion_of_total_partnership(mat_e50)
calculate_proportion_of_total_partnership(mat_e50_α)

calculate_proportion_of_individuals(mat_e50)
calculate_proportion_of_individuals(mat_e50_α)
# -

# Calculate the proportion of individuals removed from sims.
inds = [(1:11), (12:15), (16:50)]
n_ind_new = sum(mat_e50_α; dims=2) ./ KS_4W
C_ind = sum(n_ind_new)
n_inc = 0
for ind in inds[1:2]
    n_inc += (sum(mat_e50_α[ind, ind] .- mat_e50[ind, ind]; dims=2) ./ KS_4W[ind]) |> sum
end
p_red = n_inc/C_ind
println(f"Percetage of individuals removed from sims (proportion of rewired individuls): {p_red*100:.2f}%")
println(f"Remaining MSM prop: {1-p_red:.3f}")


