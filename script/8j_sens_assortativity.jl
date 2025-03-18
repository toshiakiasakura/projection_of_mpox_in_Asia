# -*- coding: utf-8 -*-
# ## Adjust assortativity based on reported Newman's assortativity coefficient

include("./utils.jl")
include("./distributions.jl")
include("./model.jl")
include("./model_intercountry.jl")

KS_4W, PK_4W = degree_and_probability_weibull(NATSAL4W_PARMS...)
scale_v = 1/4
println("Number of sexual activity group with less than 100 sexual partnership per 3 months: ", sum((KS_4W * scale_v) .< 100))

# +
inds = [(1:11), (12:15), (16:50)]
inds1 = inds[1]
inds2 = inds[2]
println(f"Check the proportion of ind1 and ind2: {sum(PK_4W[inds1]):.3f}, and {sum(PK_4W[inds2]):.3f}")

# Only monogamous individuals are modified. 
α = 0.485
α2 = 0.0

Pk_new, w1, w2 = create_rewired_degree_distribution(PK_4W, α, α2)
mat_e50_prop, mat_e50 = create_total_and_prop_partnership_dist(Pk_new, w1, w2, KS_4W, α, α2)
mat_e3 = convert_50comp_to_EricChow2016_bins(mat_e50, KS_4W)
Newman_assortativity_coef(mat_e3)
# -

# checking the original degree distribution and rewired degree distribution. 
pl = plot()
plot!(pl, 1:50, PK_4W, label="Pk_original")
plot!(pl, 1:50, Pk_new, label="Pk_new")
plot!(pl, 1:50, (1-α).*Pk_new, label="Pk_new*(1-α)")

clim = (0, maximum(mat_e50))
pl1 = heatmap(mat_e50_prop, 
    title="Proportionate mixing layer", xlabel="Sexual acitvity level", ylabel="Sexual activity level", 
    clim=clim)
pl2 = heatmap(mat_e50, title="Total layer", xlabel="Sexual acitvity level", clim=clim)
plot(pl1, pl2, size=(1000,400), dpi=300, bottom_margin=5Plots.mm, left_margin=5Plots.mm)

# Check the maximum index -> max_ind_k = 42 (not changed from before). 
println("MSM pop is reduced to $(1-α)")
targetdata, N0_MSM = read_constant_values()
Ns = N0_MSM .* (1- α) .* Pk_new
Ns[1:43]


