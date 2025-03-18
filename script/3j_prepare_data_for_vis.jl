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

include("./utils.jl")
include("./vis_util.jl")
include("./model.jl")
include("./model_intercountry.jl")

# # Data preparation
# The repository does not contain intermediate files or result files from simulations due to a large file size. 
# Fitting and simulations (which take hours to days) should be run beforehand to proceed with this file. 

# ## Filter the file

# To reduce the computation burden later on, move eligible files to another directory.
move_eligible_file_to_directory(VARS_SC1, VARS_SC1_FIL)

record_summary_statistics(VARS_SC1)

# With the same conditions as the main analysis (fil_week = 24). 
move_eligible_file_to_directory(VARS_SC1_D0415, VARS_SC1_D0415_FIL)

# With the same conditions as the main analysis (fil_week = 24). 
record_summary_statistics(VARS_SC1_D0415)

# ## Save Inc and imp. prob.

save_inc_imp(VARS_SC1_FIL)
save_inc_imp(VARS_SC2_FIL)
save_inc_imp(VARS_SC3_FIL)
save_inc_imp(VARS_SC4_FIL)

save_filtered_imp(VARS_SC1_FIL)
save_filtered_imp(VARS_SC2_FIL)
save_filtered_imp(VARS_SC3_FIL)
save_filtered_imp(VARS_SC4_FIL)

# For other sensitivity analyses.
save_inc_imp(VARS_SC1_CUT1000_FIL)
save_inc_imp(VARS_SC1_D0415_FIL)
save_inc_imp(VARS_SC1_Y2023_FIL)
save_inc_imp(VARS_SC1_ASSORT_FIL)

# For predict 1st importation dates.
save_inc_imp(VARS_SC1_10times_FIL)
save_inc_imp(VARS_SC2_10times_FIL)
save_inc_imp(VARS_SC3_10times_FIL)
save_inc_imp(VARS_SC4_10times_FIL)

# ## Generation based Sankey diagram, 0th (Japan), 1st, 2nd, 3rd.

df_upd = gen_based_Sankey_diagram(VARS_SC1_FIL, INTER_SIM_BASE.country_dict)
nothing

# ## Conditional params
# Note: This section may cause an error if individual simulation data are lacking.

cut_off = 10
tp = "Korea"

conditional_beta_SAR(VARS_SC1_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC2_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC3_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC4_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC1_CUT1000_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC1_D0415_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC1_Y2023_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC1_Y2023_FIL, cut_off, tp)

conditional_beta_SAR(VARS_SC1_ASSORT_FIL, cut_off, tp)

# # Visualisation of results

# ## Sensitivity analysis: different sexual distribution and infectious period.

var_lis =  [VARS_SC1_FIL,  VARS_SC2_FIL,  VARS_SC3_FIL,  VARS_SC4_FIL]
load_I_inc!.(var_lis)
nothing

title1 = "Natsal 4-week \n IP of 10 days"
title2 = "Natsal 1-year \n IP of 10 days"
title3 = "Natsal 4-week \n IP of 7 days"
title4 = "Natsal 4-week \n IP of 14 days"
titles_natsal = [title1, title2, title3, title4]
nothing

include("./vis_util.jl")
plot_imp_prob_global_fs_num_cnts(var_lis, titles_natsal; kwds=BAR_KWDS)

# ## Different conditionality

print_number_of_valid_sims(VARS_SC1)

cut_off = 10
title1 = "Unconditional\n"
title2 = "≥$(cut_off) cases \nin Republic of Korea"
title3 = "≥$(cut_off) cases \nin Taiwan"
title4 = "≥$(cut_off) cases \nin China"
titles_cond = [title1, title2, title3, title4]
df_mer = plot_imp_prob_global_fs_num_cnts_large_data(
    VARS_SC1, titles_cond;
    sort_col=:prob1, kwds=BAR_KWDS)

# ## Importation date

var_lis_10times =  [VARS_SC1_10times_FIL,  VARS_SC2_10times_FIL,  VARS_SC3_10times_FIL,  VARS_SC4_10times_FIL]
load_I_inc!.(var_lis_10times)

plot_1st_importation_scenarios(var_lis_10times, titles_natsal)

# ## Incidence curve

plot_multiple_incidence_curve(VARS_SC1)

# # Sensitivity analyses

# ## Sensitivity analysis: Limit maximum k to be less than 1000

plot(VARS_SC1_CUT1000_FIL)

# ## Sensitivity analysis: behavioural change

include("vis_util.jl")

n_week = SCENARIO1_D0415.mixin["targetdata"] |> length
plot_multiple_incidence_curve(VARS_SC1_D0415; 
    ylim=[0.9, 6000], anno_pos_x=5, anno_pos_y=2000, exclude_obs_panel=true, vspan=[-10, n_week])

plot(VARS_SC1_D0415_FIL; legend=(0.7, 0.1))

# ## Sensitivity analysis: international flight volume in 2023

plot(VARS_SC1_Y2023_FIL)

# ## Sensitivity analysis: increased assortativity

plot(VARS_SC1_ASSORT_FIL)
