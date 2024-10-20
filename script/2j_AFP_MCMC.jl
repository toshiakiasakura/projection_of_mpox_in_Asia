# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Julia 1.8.3
#     language: julia
#     name: julia-1.8
# ---

include("./util.jl")
include("./distributions.jl")
include("./model.jl")
include("./model_intercountry.jl")
Threads.nthreads()
theme(:ggplot2)

iterations = 12_000
n_sim_samples = 10
nparticles = 300
n_burn=2000

sim_jpn_1y_IP10 = "../tmp_afp/natsal1y_inf10_sim_jpn.jld2"
sim_jpn_1y_IP7 = "../tmp_afp/natsal1y_inf7_sim_jpn.jld2"
sim_jpn_1y_IP14 = "../tmp_afp/natsal1y_inf14_sim_jpn.jld2"
sim_jpn_4w_IP10 = "../tmp_afp/natsal4w_inf10_sim_jpn.jld2"
sim_jpn_4w_IP7 = "../tmp_afp/natsal4w_inf7_sim_jpn.jld2"
sim_jpn_4w_IP14 = "../tmp_afp/natsal4w_inf14_sim_jpn.jld2"

# # Fit and obtain japan's trajectories

# Natsal 1y data, 10 day IP. 
path = fit_β(;
    α=NATSAL1Y_PARMS[1], κ=NATSAL1Y_PARMS[2], γ2=1/10, 
    days=8, max_ind_k=50,
    iterations=iterations,
    file_prefix="natsal1y_posterior_inf7",
    σ_β=0.1,
)

# Natsal 4w data, 10 day IP. 
path = fit_β(;
    α=NATSAL4W_PARMS[1], κ=NATSAL4W_PARMS[2], γ2=1/10, 
    days=8, max_ind_k=42,
    iterations=iterations,
    file_prefix="natsal4w_posterior_inf10",
    σ_β=0.3,
)

# Natsal 4w data, 7 day IP. 
path = fit_β(;
    α=NATSAL4W_PARMS[1], κ=NATSAL4W_PARMS[2], γ2=1/7, 
    days=8, max_ind_k=42,
    iterations=iterations,
    file_prefix="natsal4w_posterior_inf10",
    σ_β=0.3,
)

# Natsal 4w data, 14 day IP. 
path = fit_β(;
    α=NATSAL4W_PARMS[1], κ=NATSAL4W_PARMS[2], γ2=1/14, 
    days=8, max_ind_k=42,
    iterations=iterations,
    file_prefix="natsal4w_posterior_inf10",
    σ_β=0.3,
)

# ## Obtain Jpana's trajectories

path1_IP10_trace = "../tmp_afp/natsal1y_posterior_inf10_20240127_032905.jls"
path5_IP10_trace = "../tmp_afp/natsal4w_posterior_inf10_20240127_021341.jls"
path5_IP7_trace = "../tmp_afp/natsal4w_posterior_inf10_20240303_140221.jls"
path5_IP14_trace = "../tmp_afp/natsal4w_posterior_inf10_20240303_133358.jls"

sim_jpn_dir_4w_IP10 = "../tmp_afp/natsal4w_inf10_sim_jpn"
sim_jpn_dir_4w_IP7 = "../tmp_afp/natsal4w_inf7_sim_jpn"
sim_jpn_dir_4w_IP14 = "../tmp_afp/natsal4w_inf14_sim_jpn"
sim_jpn_dir_1y_IP10 = "../tmp_afp/natsal1y_inf10_sim_jpn"

kwds = (n_sim_samples=n_sim_samples, nparticles=nparticles, n_burn=n_burn)

path = path5_IP10_trace
obtain_jpn_trajectories_for_each_β(
    path, sim_jpn_dir_4w_IP10; kwds...
)

path = path1_IP10_trace
obtain_jpn_trajectories_for_each_β(
    path, sim_jpn_dir_1y_IP10; kwds...
)

path = path5_IP7_trace
obtain_jpn_trajectories_for_each_β(
    path, sim_jpn_dir_4w_IP7; kwds...
)

path = path5_IP14_trace
obtain_jpn_trajectories_for_each_β(
    path, sim_jpn_dir_4w_IP14; kwds...
)

# ### Run inter-country simulations 

n_inter_sim = 50_000

include("model_intercountry.jl")

run_and_save_intercountry_model(sim_jpn_dir_4w_IP10; Korea_cond=false)

run_and_save_intercountry_model(sim_jpn_dir_1y_IP10; Korea_cond=true)

run_and_save_intercountry_model(sim_jpn_dir_4w_IP7; Korea_cond=true)

run_and_save_intercountry_model(sim_jpn_dir_4w_IP14; Korea_cond=true)



# ## Thinned samples

chn = show_traceplot(path5_IP10_trace)

chn = show_traceplot(path1_IP10_trace)

chn = show_traceplot(path5_IP7_trace)

chn = show_traceplot(path5_IP14_trace)


