# -*- coding: utf-8 -*-
include("./utils.jl")
include("./utils_AIMSS.jl")
include("./distributions.jl")


# ## Cleaning process

path = "../data_AIMSS/AIMSS.csv" 
df = CSV.read(path, DataFrame)
# Rename countries to match the names in the manuscript. 
rep_cnt = Dict("China" => "China (mainland)",  "Korea (South)" => "Republic of Korea")
df = @transform(df, :countryname = replace.(:countryname, rep_cnt...))
size(df)

df_fil2 = clean_process1(df)
df_tmp = rep_and_parse_values(df_fil2)
# Here, we do not accept contradicted answers in the dataset.
df_ana = @subset df_tmp begin
    .!((:nmsex_num .< :nregular_num) )
    .!((:nmsex_num .< :ncasual_num) )
    .!((:nmsex_num .< :ncomm_num) )
end
size(df_ana) |> println
df_ana[:, :countryname] |> countmap

println("Number of countries reported in this study: ", df_ana[:, :countryname] |> unique |> length)

# # Visualisation 

include("./utils_AIMSS.jl")
tab_long1 = prepare_tab_long1(df_ana)
tab_long2 = prepare_tab_long2(df_ana)
nothing

# +
color = reshape(
    [palette(:ice)[end:-1:begin][31:42:(12*21)]..., 
        palette(:algae)[36:42:(11*21)]...,
        :red1, :blue, :cyan,
    ], 
    (1, 14)
)
kwds_bar = (
        foreground_color_legend = nothing,
        background_color_legend = nothing,
        ylabel = "Proportion of each category (%)",
        xlabel = "Number of sexual partners in the previous 6 months",
        ylim=[0, 80],
        color=color,
        markersize=2.5,
)

pl1 = plot_AIMSS_bar1(tab_long1, kwds_bar)
pl2 = plot_AIMSS_bar2(tab_long2, kwds_bar)
annotate!(pl1, (0.05, 0.93), text("A", :black, :left, :bottom, 32))
annotate!(pl2, (0.05, 0.93), text("B", :black, :left, :bottom, 32))
plot(pl1, pl2, layout=(2,1), dpi=300, size=(800, 1000),
    left_margin=5Plots.mm, top_margin=5Plots.mm
)
# -










