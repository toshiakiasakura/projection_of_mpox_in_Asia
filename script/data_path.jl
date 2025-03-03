# Data/results path summaries.

PATH_LINE_MASTER = "../data/JPN_linelist_master_20230707.csv"
PATH_NIID = "../data/jpn_epicurve_NIID.csv"
PATH_FLIGHT = "../data/flight/selected_flight_matrix.csv"
PATH_POP = "../data/pop_size_edit.csv"
PATH_ISO_CODE = "../data/iso2_iso3.csv"
PATH_OWID = "../data/owid-mpox-data-20250121.csv"

PATH_JPN_0707 = "../tmp_fix_results/mpox_japan_weekly_data_0707.csv"
PATH_JPN_0415 = "../tmp_fix_results/mpox_japan_weekly_data_0415.csv"

# Path for estimated beta calibrated to the Japanese incidence data. .
PATH_SC1 = "../tmp_afp/natsal4w_posterior_inf10_20241021_053947.jls" # 12_000
PATH_SC2 = "../tmp_afp/natsal1y_posterior_inf10_20241021_095129.jls"
PATH_SC3 = "../tmp_afp/natsal4w_posterior_inf7_20241021_073031.jls"
PATH_SC4 = "../tmp_afp/natsal4w_posterior_inf14_20241021_070357.jls"
PATH_SC1_CUT1000 = "../tmp_afp/natsal4w_posterior_inf10_cut1000_20241102_155228.jls"
PATH_SC1_D0415 = "../tmp_afp/natsal4w_posterior_inf10_d0415_20241107_100824.jls"
PATH_SC1_ASSORT = "../tmp_afp/natsal4w_posterior_inf10_assort_20250228_143847.jls"

# Path for trajectories in Japan, from the estimated beta values.
SIM_JPN_DIR_SC1 = "../tmp_afp/natsal4w_inf10_sim_jpn"
SIM_JPN_DIR_SC2 = "../tmp_afp/natsal1y_inf10_sim_jpn"
SIM_JPN_DIR_SC3 = "../tmp_afp/natsal4w_inf7_sim_jpn"
SIM_JPN_DIR_SC4 = "../tmp_afp/natsal4w_inf14_sim_jpn"
SIM_JPN_DIR_SC1_CUT1000 = "../tmp_afp/natsal4w_inf14_sim_jpn_cut1000"
SIM_JPN_DIR_SC1_D0415= "../tmp_afp/natsal4w_inf14_sim_jpn_d0415"
SIM_JPN_DIR_SC1_ASSORT = "../tmp_afp/natsal4w_inf10_sim_jpn_assort"

# DIR for intercountry spread. # 12_000 .
DIR_INTERCNT_SC1 = "../tmp_results/20241103_105517_sc1"
DIR_INTERCNT_SC1_FIL = "$(DIR_INTERCNT_SC1)_fil"
DIR_INTERCNT_SC2 = "../tmp_results/20241103_101351_sc2"
DIR_INTERCNT_SC3 = "../tmp_results/20241103_145415_sc3"
DIR_INTERCNT_SC4 = "../tmp_results/20241103_153744_sc4"
DIR_INTERCNT_SC1_CUT1000 = "../tmp_results/20241103_162220_sc1_cut1000"
DIR_INTERCNT_SC1_D0415 = "../tmp_results/20241122_000749_sc1_d0415"
DIR_INTERCNT_SC1_D0415_FIL = "$(DIR_INTERCNT_SC1_D0415)_fil"
DIR_INTERCNT_SC1_Y2023_FIL = "../tmp_results/20241112_004906_sc1_y2023"
DIR_INTERCNT_SC1_ASSORT_FIL = "../tmp_results/20250228_182302_sc1_assort_fil"

DIR_INTERCNT_SC1_10times_FIL = "../tmp_results/20250218_000032_sc1_10times"
DIR_INTERCNT_SC2_10times_FIL = "../tmp_results/20250218_060009_sc2_10times"
DIR_INTERCNT_SC3_10times_FIL = "../tmp_results/20250209_134814_sc3_10times"
DIR_INTERCNT_SC4_10times_FIL = "../tmp_results/20250218_122406_sc4_10times"

# Path and directory for summarising results.
Base.@kwdef mutable struct VisVars
	path_sim_res::String
	suffix::String
	path_inc::String = "../tmp_results/I_inc_$(suffix).jls"
	path_imp::String = "../tmp_results/imp_prob_$(suffix).csv"
	path_imp_fil::String = "../tmp_results/imp_prob_$(suffix)_fil.csv"
	path_imp_date::String = "../tmp_results/imp_dates_$(suffix).jls"
	path_imp_date_fil::String = "../tmp_results/imp_dates_$(suffix)_fil.jls"
	path_exp_imp_gen::String = "../tmp_results/exp_imp_gen_pair_$(suffix).csv"
	path_exp_imp_gen_fil::String = "../tmp_results/exp_imp_gen_pair_$(suffix)_fil.csv"
    path_rec::String = "../tmp_results/summary_quantity_$(suffix).jld2"
	I_inc::Array{Float64, 3} = fill(Inf, 1, 1, 1)
end

# Before filtering
VARS_SC1 = VisVars(path_sim_res = DIR_INTERCNT_SC1, suffix = "path_sc1")
VARS_SC1_D0415 = VisVars(path_sim_res = DIR_INTERCNT_SC1_D0415, suffix = "path_sc1_d0415")

# FIL represents the condition on Korea incidence.
VARS_SC1_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC1_FIL, suffix = "path_sc1_fil")
VARS_SC2_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC2, suffix = "path_sc2_fil")
VARS_SC3_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC3, suffix = "path_sc3_fil")
VARS_SC4_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC4, suffix = "path_sc4_fil")
VARS_SC1_CUT1000_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC1_CUT1000,
	suffix = "path_sc1_cut1000_fil")
VARS_SC1_D0415_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC1_D0415_FIL,
    suffix = "path_sc1_d0415_fil")
VARS_SC1_Y2023_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC1_Y2023_FIL,
    suffix = "path_sc1_y2023_fil")
VARS_SC1_ASSORT_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC1_ASSORT_FIL,
	suffix = "path_sc1_assort_fil")

VARS_SC1_10times_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC1_10times_FIL,
	suffix = "path_sc1_10times_fil")
VARS_SC2_10times_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC2_10times_FIL,
	suffix = "path_sc2_10times_fil")
VARS_SC3_10times_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC3_10times_FIL,
	suffix = "path_sc3_10times_fil")
VARS_SC4_10times_FIL = VisVars(path_sim_res = DIR_INTERCNT_SC4_10times_FIL,
	suffix = "path_sc4_10times_fil")
# Note: This file is not uploaded to the repository.
# TODO: consider the removal of this path or not.
PATH_OAG_RATIO = "../data/OAG_summary/OAG_ratio_2023to2019_5th_95th_truncated.jld2"

