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

include("./utils.jl")
include("./distributions.jl")
include("./model.jl")
include("./figure_utils.jl")
include("./data_path.jl")

# # Explore Japan mpox cases

# +
df = CSV.read(PATH_LINE_MASTER, DataFrame)
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

haskey(Dict(:a=>"b"), :a)

# +
# For saving purpose.
df_fil = @subset df $med_date .>= Date(2023, 1, 1)
tab_7d = convert_to_weekly_case(df_fil, med_date)
CSV.write(PATH_JPN_0707, tab_7d)
println("Number of week until 0707: ", nrow(tab_7d))

df_fil_0415 = @subset df_fil $med_date .<= Date(2023, 4, 15)
tab_7d_0415 = convert_to_weekly_case(df_fil_0415, med_date)
CSV.write(PATH_JPN_0415, tab_7d_0415)
println("Number of week until 0415: ", nrow(tab_7d_0415))
# -

Japan_weekly_epicurve_2022_2023(df)

# # Sexual partner distribution.

visualise_sexual_partner_distribution(KS_4W, PK_4W, KS_1Y, PK_1Y)

visualise_sexual_partner_distribution(KS_4W, PK_4W_CUT1000, KS_1Y, PK_1Y)
