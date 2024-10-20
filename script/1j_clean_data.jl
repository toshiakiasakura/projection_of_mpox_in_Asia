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

# # Fitting to Japanese cases

using Pkg
Pkg.instantiate()

include("./distributions.jl")
include("./model.jl")
include("./util.jl")
include("./figure_utils.jl")

# # Explore Japan mpox cases

# +
path = "../data/JPN_linelist_master_20230707.csv"
df = CSV.read(path, DataFrame)
names(df) |> println
med_date_pre = "Date of medical attendance"
med_date = "medical_attendance"
date_repo = "Date of reporting from MHLW"

df = @chain df begin
    @transform(:count = 1,
        $date_repo = Date.($date_repo, "yyyy/mm/dd"),
    )
end
df[!, med_date] = ifelse.( ismissing.(df[:, med_date_pre]),
    df[:, date_repo] .- Dates.Day(9),
    passmissing(Date).(df[:, med_date_pre], "yyyy/mm/dd")
    )
nothing
# -

# For saving purpose.
df_fil = @subset df $med_date .>= Date(2023, 1, 1) 
tab_7d = convert_to_weekly_case(df_fil, med_date)
CSV.write("../tmp_fix_results/mpox_japan_weekly_data_0707.csv", tab_7d)
Japan_weekly_epicurve_2022_2023()

# # Sexual partner distribution.

# TODO: For the review purpose, some cut-off values should be set and adjusted the degree and probability. 
ks_1y, Pk_1y = degree_and_probability_weibull(NATSAL1Y_PARMS...) # Natsal 1 year.
ks_4w, Pk_4w = degree_and_probability_weibull(NATSAL4W_PARMS...) # Natsal 4 weeks.
nothing

visualise_sexual_partner_distribution(ks_1y, Pk_1y, ks_4w, Pk_4w)



