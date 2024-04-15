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

# # Explore Japan mpox cases

path = "../data/JPN_linelist_master_20230707.csv"
df = CSV.read(path, DataFrame)
names(df) |> println
med_date_pre = "Date of medical attendance"
med_date = "medical_attendance"
date_repo = "Date of reporting from MHLW"

df[!, date_repo] .= Date.(df[:, date_repo], "yyyy/mm/dd")
df[:, date_repo] .- Dates.Day(9)
df[!, med_date] = ifelse.( ismissing.(df[:, med_date_pre]),
    df[:, date_repo] .- Dates.Day(9),
    passmissing(Date).(df[:, med_date_pre], "yyyy/mm/dd")
    )
nothing

# +
# (df[:, date_repo] .- df[:, med_date])  .|> Dates.days |> mean
# 8.66120218579235
# -

df[:, :count] .= 1
tab = @pipe groupby(df, med_date) |> combine(_, :count => sum)
bar(tab[:, med_date], tab[:, :count_sum], size=(900,400), label=:none,
    xlabel="Medical attendance date",
    ylabel="Daily incidence",
    fmt=:png,
    left_margin=5Plots.mm,
    bottom_margin=5Plots.mm
)

cond1 = df[:, med_date] .>= Date(2023,1,1)
#cond2 = df[:, "Prefecture reported"] .!= "Osaka"
#dfM = df[cond1 .& cond2, :]
dfM = df[cond1, :]
tab = @pipe groupby(dfM, med_date) |> combine(_, :count => sum)
tab = sort(tab, [med_date])
pl = bar(tab[:, med_date], tab[:, :count_sum], size=(900,400),
    #title="Check epicurve for mpox in Japan.",
    label=nothing,
    dpi=300,
    fmt=:png,
    xtickfontsize=10,
    ytickfontsize=10,
    guidefontsize=14,
    xlabel="Medical attendance date",
    ylabel="Daily count",
    left_margin=15Plots.pt,
    bottom_margin=15Plots.pt,
    foreground_color_legend=nothing,
)
#savefig(pl, "../fig/epi_curve_JPN.png")
display(pl)

# ### Data preparation

"""Given dataframe, it tabulates to calculate the weekly numbers.
"""
function convert_to_weekly_case(df_::DataFrame, col::String)
    tab = @pipe groupby(df_, col) |> combine(_, :count => sum)
    sort!(tab, med_date)
    tab_7d = weekly_num_cases(tab, col)
    return tab_7d
end

# +
df[:, :count] .= 1

# For saving purpose.
cond1 = (df[:, med_date] .>= Date(2023,1,1))
tab_7d = convert_to_weekly_case(df[cond1, :], med_date)
CSV.write("../tmp_fix_results/mpox_japan_weekly_data_0707.csv", tab_7d)
# +
# Obtain from 2022-07-25
tab_7d = convert_to_weekly_case(df, med_date)

cond = coalesce.(df[:, "Travel history"], "None") .!= "None"
tab_7d_travel = convert_to_weekly_case(df[cond, :], med_date)
nothing
# -

tab_vis = leftjoin(tab_7d, tab_7d_travel, on=:timestamp, makeunique=true)
tab_vis[:, :count_sum_1] = coalesce.(tab_vis[:,:count_sum_1], 0)
tab_vis[:, :no_history] = tab_vis[:, :count_sum] .- tab_vis[:, :count_sum_1]
date = tab_vis[:, :timestamp]
cnts_no_history = tab_vis[:, :count_sum]
cnts_history = tab_vis[:, :count_sum_1]
nothing

xtick_label = [Date(2022,7,25), Date(2022, 10, 10), Date(2023,1,2), Date(2023,3,27), Date(2023,6,19)]
pl = bar(date, cnts_no_history, label="No travel history")
bar!(pl, date, cnts_history, label="Travel history")
plot!(pl,
    label=:none,
    dpi=300,
    fmt=:png,
    xtickfontsize=10,
    ytickfontsize=10,
    guidefontsize=14,
    xlabel="Medical attendance date",
    ylabel="Weekly incidence",
    foreground_color_legend=nothing,
    background_color_legend=nothing,
    left_margin=15Plots.pt,
    bottom_margin=15Plots.pt,
    size=(800,400),
    xticks = (xtick_label, xtick_label),
    ylim=[0,20]
)
savefig(pl, "../fig/epi_curve_JPN.png")
#display(pl)

# # Sexual partner distribution.

ks_1y, Pk_1y = degree_and_probability_weibull(0.10, 0.77) # Natsal 1 year.
ks_4w, Pk_4w = degree_and_probability_weibull(0.16, 0.87529) # Natsal 4 weeks.
nothing

# +
x = 1:50
rate_kwds = (yaxis=:log10, label=:none,
    xlim=[0.5,51], ylim=[1e-3, 35], 
    marker=(:circle, 3),
    # fillrange=1e-3,
    yticks = ([1e-2, 1e-1, 1.0, 10], ["0.01", "0.1", "1.0", "10"]),
    bar_edges=false, bar_width=1.0, lw=0.3,
    xlabel="Sexual activity group",
)
prob_bar_kwds = (yaxis=:log10, label=:none, ylim=[1e-9, 0.5],
    fillrange=1e-9,
    bar_width=1.0, lw=0.3,
    ylabel="Proportion of \neach sexual activity group",
    xlabel="Sexual activity group",
    xlim=[0.5,50.5],
)

pl1 = plot(x, ks_4w/365;
    ylabel="Rate of sexual encounters \nper day",
    rate_kwds...)
annotate!(pl1, 5, 10, text("A", :black, :left, 18))

pl2 = bar(x, Pk_4w; #title="Natsal 4 week",
    prob_bar_kwds...)
annotate!(pl2, 40, 8*1e-2, text("B", :black, :left, 18))

pl3 = bar(x, Pk_1y; #title="Natsal 1 year",
    prob_bar_kwds...)
annotate!(pl3, 40, 8*1e-2, text("C", :black, :left, 18))
#pl4 = bar(x, ks_1y/365; rate_kwds...)
layout = @layout [a b c]
pl = plot(pl1, pl2, pl3, layout=layout,
    dpi=300, size=(1000,300),
    left_margin=10Plots.mm,
    bottom_margin=7Plots.mm,
)
savefig(pl, "../fig/sexual_distribution.png")
# -




