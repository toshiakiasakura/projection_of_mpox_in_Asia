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
include("./vis_util.jl")
include("./model.jl")
include("./model_intercountry.jl")

# # Data preparation
# Note: Files in `../tmp_results/**` is created in 2j, and those file size is around 900GB for Natsal 4w IP10 (leaving all simulation resutls), and around 15GB for other scenarios (leaving only conditioned results).

path_flight = "../data/flight/selected_flight_matrix.csv"
path_pop = "../data/pop_size_edit.csv"
m_return, N0_MSM, country_dict = read_inter_country_data(path_flight, path_pop)

# +
vars1_IP10 = VisVars(
    path_sim_res = "../tmp_results/20240411_003226", # 50_000, conditional only with β,
    suffix = "path1_IP10"
)

vars5_IP10 = VisVars(
    path_sim_res = "../tmp_results/20240411_005101_fil", # 50_000 with Korea condition.
    suffix = "path5_IP10"
)

vars5_IP7 = VisVars(
    path_sim_res = "../tmp_results/20240411_013445", # 50_000, only conditional one with β,
    suffix = "path5_IP7"
)

vars5_IP14 = VisVars(
    path_sim_res = "../tmp_results/20240411_005430", # 50_000, only conditional one with β,
    suffix = "path5_IP14"
)
# -

# ## Filter the file

vars5_IP10_ori = VisVars(
    path_sim_res = "../tmp_results/20240411_005101", # 50_000 without any conditions.
    suffix = "path5_IP10"
)
path_rec = "../tmp_results/summary_quantity.jld2"

move_eligible_file_to_directory(vars5_IP10_ori, "../tmp_results/20240411_005101_fil")

record_summary_statistics(vars5_IP10_ori, path_rec)

record_summary_statistics(vars5_IP10_ori, "../tmp_results/summary_quantity_test.jld2"; nmax=50000)

# # Save Inc and imp. prob.

nmax=50_000

save_inc_imp(vars5_IP10)

save_inc_imp(vars1_IP10)
save_inc_imp(vars5_IP7)
save_inc_imp(vars5_IP14)
save_inc_imp(vars5_IP10)

save_filtered_imp(vars1_IP10)
save_filtered_imp(vars5_IP7)
save_filtered_imp(vars5_IP14)
save_filtered_imp(vars5_IP10)

# ## Generation based Sankey diagram, 0th (Japan), 1st, 2nd, 3rd.

df_upd = gen_based_Sankey_diagram(vars5_IP10.path_sim_res)
CSV.write(vars5_IP10.path_exp_imp_gen, df_upd)

df_fil = filter_exp_gen(vars5_IP10)
CSV.write(vars5_IP10.path_exp_imp_gen_fil, df_fil)

# # Conditional params
# Note: This section may cause an error if individual simulation data are lacking.

include("model_intercountry.jl")
include("vis_util.jl")

cut_off = 10
tp = "Korea"

conditional_beta_SAR(vars5_IP10, cut_off, tp; n_sim=50_000)

conditional_beta_SAR(vars1_IP10, cut_off, tp; n_sim=50_000)

conditional_beta_SAR(vars5_IP7, cut_off, tp; n_sim=50_000)

conditional_beta_SAR(vars5_IP14, cut_off, tp; n_sim=50_000)
