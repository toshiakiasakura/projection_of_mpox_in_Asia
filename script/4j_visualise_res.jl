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
# The files needed to visualise figures are uploaded. Skip the following `Filter the file`, `Save Inc and imp. prob`, and `Generation based Sankey diagram, 0th (Japan), 1st, 2nd, 3rd.` sections. 
#
# Note: Files in `../tmp_results/**` is created in 2j, and those file size is around 900GB for Natsal 4w IP10 (leaving all simulation resutls), and around 15GB for other scenarios (leaving only conditioned results).

β = 0.091
SAR = 1- exp(-β)

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
path_rec = "../tmp_results/summary_quantity.jld2"
# -

# # Visualisation

include("vis_util.jl")
ind0_cnt = 15

load_I_inc!(vars1_IP10)
load_I_inc!(vars5_IP7)
load_I_inc!(vars5_IP14)
load_I_inc!(vars5_IP10)
nothing

# ## Figure preparation

labelfontsize = 8

# Porb. importation risks.
n_cnt = 41
yticks_empty = (1:n_cnt, ["" for i in 1:n_cnt])
kwds = Dict(
    :xlim => [0, 105],
    :orientation=>:horizontal,
    :ylim => [0, n_cnt+1],
    :label=>"",
    :left_margin=>7Plots.pt,
    :right_margin=>0Plots.pt,
    :upper_margin=>40Plots.pt,
    :xlabel=>"Simulated importation\nprobability (%)",
    :titlefontsize=>10,
    :ytickfontsize=> labelfontsize,
    :xlabelfontsize=>labelfontsize,
    :ylabelfontsize=>labelfontsize,
)

title1 = "Natsal 4-week \n IP of 10 days"
title2 = "Natsal 1-year \n IP of 10 days"
title3 = "Natsal 4-week \n IP of 7 days"
title4 = "Natsal 4-week \n IP of 14 days"
visualise_imp_prob_global_fs_num_cnts(
    vars5_IP10.I_inc,
    vars1_IP10.I_inc,
    vars5_IP7.I_inc,
    vars5_IP14.I_inc,
    title1, title2, title3, title4;
    kwds=kwds
)

# # Different conditionality

include("model_intercountry.jl")
include("vis_util.jl")

#path_test = "../tmp_results/summary_quantity_test.jld2"
# path_rec = "../tmp_results/summary_quantity.jld2"
sum_obj = load(path_rec)
r_obj = sum_obj["r_obj"]
nothing

println("# of all sims: ", length(r_obj.fs))
println("# of Korea: ", sum(sum_obj["Korea_ind"]))
println("# of Taiwan: ", sum(sum_obj["Taiwan_ind"]))
println("# of China: ", sum(sum_obj["China_ind"]))

cut_off = 10
title1 = "Unconditional\n"
title2 = "≥$(cut_off) cases \nin South Korea"
title3 = "≥$(cut_off) cases \nin Taiwan"
title4 = "≥$(cut_off) cases \nin China"
df_mer = visualise_imp_prob_global_fs_num_cnts_large_data(
    path_rec, title1, title2, title3, title4;
    sort_col=:prob1, kwds=kwds)

# # Importation date

include("vis_util.jl")

# +
I_inc1 = @pipe vars5_IP10.I_inc |> filter_I_inc(_; tp="Korea", cut_off=cut_off)
I_inc2 = @pipe vars1_IP10.I_inc |> filter_I_inc(_; tp="Korea", cut_off=cut_off)
I_inc3 = @pipe vars5_IP7.I_inc  |> filter_I_inc(_; tp="Korea", cut_off=cut_off)
I_inc4 = @pipe vars5_IP14.I_inc |> filter_I_inc(_; tp="Korea", cut_off=cut_off)

df_imp1 = filtered_imp_prob(I_inc1, :prob1)
df_imp2 = filtered_imp_prob(I_inc2, :prob2)
df_imp3 = filtered_imp_prob(I_inc3, :prob3)
df_imp4 = filtered_imp_prob(I_inc4, :prob4)

df_mer = create_imp_prob_dataframe(df_imp1, df_imp2, df_imp3, df_imp4;
    sort_col=:prob1)
df_mer = insert_obs_date(df_mer)
sort!(df_mer, :prob1)
nothing
# -

imp_p5_IP10 = quantile_importation_dates(I_inc1)
imp_p1_IP10 = quantile_importation_dates(I_inc2)
imp_p5_IP7 = quantile_importation_dates(I_inc3)
imp_p5_IP14 = quantile_importation_dates(I_inc4)
nothing

pl5_IP10 = imp_date_vis(imp_p5_IP10, df_mer; title="Natsal 4-week, IP of 10 days", yticks=true)
pl1_IP10 = imp_date_vis(imp_p1_IP10, df_mer; title="Natsal 1-year, IP of 10 days", yticks=false, label=false)
plot!(pl5_IP10, left_margin=15Plots.mm, right_margin=2Plots.mm)
plot!(pl1_IP10, left_margin=0Plots.mm, right_margin=5Plots.mm)
pl5_IP7 = imp_date_vis(imp_p5_IP7, df_mer  ; title="Natsal 4-week, IP of 7 days", yticks=true, label=false)
pl5_IP14 = imp_date_vis(imp_p5_IP14, df_mer; title="Natsal 4-week, IP of 14 days", yticks=false, label=false)
#plot!(pl5_IP14, left_margin=-20Plots.mm)
nothing

plot([pl5_IP10, pl1_IP10, pl5_IP7, pl5_IP14]...,
    layout = @layout[a b; c d] , size=(1000, 1500),
    fmt=:png,
)

# ## Incidence curve

function weekly_incidence_figure(jpn_weekly::Matrix, title;
        xlabel="", ylabel="",
        ylim=[0.9, 30], yticks
    )
    ps = [0.25, 0.50, 0.75, 0.95]
    #jpn_weekly = one_country_I_inc_to_weekly(I_inc[:, :, ind0_cnt])
    q_dict = quantiles_over_week(jpn_weekly, ps)
    n_week = size(jpn_weekly)[2]
    pl = plot(
        yaxis=:log10,
        lw=4,
        xlim=[0, 159],
        ylim=ylim,
        yticks=yticks,
        titlefontsize=14,
        legendfontsize=10,
        legendtitlefontsize=10,
        xtickfontsize=10,
        ytickfontsize=10,
        title=title,
        xlabel=xlabel, ylabel=ylabel,
        fmt=:png, dpi=200,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
    )
    vspan!(pl, [-10,23], color=:black, alpha=0.2, label=:none)
    for p in ps
        lw = p == 0.5 ? 2 : 1
        ls = p == 0.95 ? :dot : :solid

        q = q_dict[p]
        p_lab = Int64(p*100)
        plot!(pl, 1:n_week, q .+ 1,
            lw=lw, ls=ls, color="blue",
            label="$(p_lab)th",
            legend_title="Quantile",
        )
    end
    pl
end

# +
# Load data
sum_obj = load(path_rec)
r_obj = sum_obj["r_obj"]
# Convert data to appropriate data type.
jpn_weekly_all = mapreduce(permutedims, vcat, r_obj.jpn_weekly)
Korea_ind = [i for (i, x) in enumerate(sum_obj["Korea_ind"])
             if x == true]
jpn_weekly_Korea = jpn_weekly_all[Korea_ind, :]

# Japanese data
jpn_obs_curve = CSV.read("../data/jpn_epicurve_NIID.csv", DataFrame)
n_jpn = nrow(jpn_obs_curve)
nothing
# -

xlim=[0, 159]
ylim = [0.9, 2000]
ylabel = "Weekly cases + 1"
yticks = ([1, 10, 100, 1000], ["1", "10", "100", "1000"])
pl1 = weekly_incidence_figure(jpn_weekly_Korea, "";
    ylim=ylim, ylabel=ylabel, yticks=yticks)
pl2 = weekly_incidence_figure(jpn_weekly_all, "";
    xlabel="Epidemiological week from 16 January 2023\n(medical attendance date)",
    ylabel=ylabel,
    ylim=ylim, yticks=yticks
)
pl3 = bar(jpn_obs_curve[:,:week_onset], jpn_obs_curve[:, :case] .+ 1,
    bar_width=1.0, label="",
    xlabel="Epidemiological week from 1 January 2023 \n(symptom onset)",
    ylabel=ylabel,
    yaxis=:log10,
    xlim=xlim,
    ylim=ylim,
    xtickfontsize=10,
    ytickfontsize=10,
    yticks=yticks,
)
# "-2" is to be matched with the fitting period.
vspan!(pl3, [-10, 23 - 2], color=:black, alpha=0.2, label=:none)
annotate!(pl1, 10, 500, text("A", :black, :left, 24))
annotate!(pl2, 10, 500, text("B", :black, :left, 24))
annotate!(pl3, 10, 500, text("C", :black, :left, 24))
layout = @layout [a;b;c]
plot(pl1, pl2, pl3, layout=layout, size=(800,800), dpi=300)







