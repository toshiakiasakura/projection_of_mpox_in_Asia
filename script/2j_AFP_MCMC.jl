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

include("./data_path.jl")
include("./utils.jl")
include("./distributions.jl")
include("./model.jl")
include("./model_intercountry.jl")
theme(:ggplot2)
Threads.nthreads()

iterations = 12_000 
nparticles = 300
n_burn= 2_000 

# # Fit and obtain Japan's trajectories

kwds = (days=SIM_DAY, iterations=iterations, nparticles=nparticles) 
nothing

sc1 = fit_β(SCENARIO1; kwds...)
sc2 = fit_β(SCENARIO2; kwds...)
sc3 = fit_β(SCENARIO3; kwds...)
sc4 = fit_β(SCENARIO4; kwds...)
# 25201.544877 seconds (179.06 G allocations: 35.226 TiB, 11.89% gc time, 0.00% compilation time)
sc1_cut1000 = fit_β(SCENARIO1_CUT1000; kwds...)

kwds = (days=SIM_DAY, iterations=22_000, nparticles=nparticles) 
sc1_d0415 = fit_β(SCENARIO1_D0415; kwds...)
# 27408.667971 seconds (186.56 G allocations: 36.808 TiB, 12.79% gc time, 0.01% compilation time)
sc1_assort = fit_β(SCENARIO1_ASSORT; kwds...)

# ## Obtain Jpana's trajectories

kwds = (nparticles=nparticles, n_burn=n_burn)

# 4414.548616 seconds (19.39 G allocations: 20.963 TiB, 50.07% gc time)
@time begin obtain_jpn_trajectories_for_each_β(PATH_SC1, SIM_JPN_DIR_SC1; kwds... ) end

@time begin obtain_jpn_trajectories_for_each_β(PATH_SC2, SIM_JPN_DIR_SC2; kwds... ) end

@time begin obtain_jpn_trajectories_for_each_β(PATH_SC3, SIM_JPN_DIR_SC3; kwds... ) end

@time begin obtain_jpn_trajectories_for_each_β(PATH_SC4, SIM_JPN_DIR_SC4; kwds... ) end

@time begin obtain_jpn_trajectories_for_each_β(PATH_SC1_CUT1000, SIM_JPN_DIR_SC1_CUT1000; kwds... ) end

# 3072.047027 seconds (12.98 G allocations: 9.102 TiB, 36.97% gc time, 0.77% compilation time
@time begin obtain_jpn_trajectories_for_each_β(PATH_SC1_D0415, SIM_JPN_DIR_SC1_D0415; 
        sc=SCENARIO1_D0415, n_thin=20, kwds...) end

@time begin obtain_jpn_trajectories_for_each_β(PATH_SC1_ASSORT, SIM_JPN_DIR_SC1_ASSORT; kwds... ) end

# ### Run inter-country simulations 

include("model_intercountry.jl")
# NOTE: This code are not executable since the flight data is not uploaded to the repository.
INTER_SIM_BASE = return_inter_sim_base()

@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC1, INTER_SIM_BASE;  Korea_cond=false) end 

@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC1, INTER_SIM_BASE;  Korea_cond=true) end 

@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC2, INTER_SIM_BASE; Korea_cond=true) end 

@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC3, INTER_SIM_BASE; Korea_cond=true) end 

@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC4, INTER_SIM_BASE; Korea_cond=true) end 

@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC1_CUT1000, INTER_SIM_BASE; Korea_cond=true) end 

@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC1_D0415, INTER_SIM_BASE; Korea_cond=false) end 

# TODO: Check OAG results to be left here or not.
INTER_SIM_2023 = return_inter_sim_base()
r_mul = JLD2.load(PATH_OAG_RATIO)["r_mul"]
INTER_SIM_2023.m_return = INTER_SIM_2023.m_return .* r_mul
# 2450.774791 seconds (2.67 G allocations: 1.615 TiB, 15.35% gc time, 0.27% compilation time: 10% of which was recompilation)
@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC1, INTER_SIM_2023; Korea_cond=true) end 
nothing

INTER_SIM_2023 = return_inter_sim_base()
INTER_SIM_2023.N0_MSM = trunc.(Int64, SCENARIO1_ASSORT.mixin["red_msm_prop"] .* INTER_SIM_2023.N0_MSM)
@time begin run_and_save_intercountry_model(SIM_JPN_DIR_SC1_ASSORT, INTER_SIM_2023; Korea_cond=true) end 

# ## Thinned samples

chn = show_traceplot(PATH_SC1; n_burn=n_burn)

chn = show_traceplot(PATH_SC2; n_burn=n_burn)

chn = show_traceplot(PATH_SC3; n_burn=n_burn)

chn = show_traceplot(PATH_SC4; n_burn=n_burn)

chn = show_traceplot(PATH_SC1_CUT1000; n_burn=n_burn)

chn = show_traceplot(PATH_SC1_D0415; n_burn=n_burn, n_thin=20)

chn = show_traceplot(PATH_SC1_ASSORT; n_burn=n_burn)
