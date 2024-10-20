
Base.@kwdef mutable struct VisVars
    path_sim_res::String
    suffix::String
    path_inc::String = "../tmp_results/I_inc_$(suffix).jl"
    path_imp::String = "../tmp_results/imp_prob_$(suffix).csv"
    path_imp_fil::String = "../tmp_results/imp_prob_$(suffix)_fil.csv"
    path_imp_date::String = "../tmp_results/imp_dates_$(suffix).jl"
    path_imp_date_fil::String = "../tmp_results/imp_dates_$(suffix)_fil.jl"
    path_exp_imp_gen::String = "../tmp_results/exp_imp_gen_pair_$(suffix).csv"
    path_exp_imp_gen_fil::String = "../tmp_results/exp_imp_gen_pair_$(suffix)_fil.csv"
    I_inc::Array{Float64,3} = fill(Inf, 1, 1, 1)
end

function load_I_inc!(vars::VisVars; check=false)
    if (check == true) & (isinf(vars.I_inc[1, 1, 1]) == false)
        return nothing
    end
    vars.I_inc = JLD2.load_object(vars.path_inc)
    return nothing
end

function read_imp_prob(path, col)
    cols = [:iso_code, :location, :imp_prob]
    df_imp = @pipe CSV.read(path, DataFrame)[:, cols] |>
                   DataFrames.rename(_, :imp_prob => col) |>
                   filter(x -> x.iso_code != "JPN", _)
    df_imp
end

function create_imp_prob_dataframe(df_mer; sort_col::Symbol=:none)
    # Add Region information.
    df_UN = CSV.read("../data/UN_Asia_list.csv", DataFrame)
    df_mer = leftjoin(df_mer, df_UN[:, [:Code, :Region]], on=:iso_code => :Code)
    cond = df_mer[:, :Region] .|> ismissing
    df_mer[cond, :] |> display
    df_mer[cond, :Region] .= "Eastern Asia"
    df_mer[:, :Region] |> unique |> display

    df_mer = sort(df_mer, [sort_col], rev=false)

    # Country name adjusting
    rep_dic = Dict(
        "China" => "China (mainland)",
        "Viet Nam" => "Vietnam",
        "Iran, Islamic Republic of" => "Iran",
        "Lao People's Democratic Republic" => "Laos",
        "Brunei Darussalam" => "Brunei",
        "Syrian Arab Republic" => "Syria",
        "Hong Kong, China" => "Hong Kong",
        "Korea, Republic of" => "Republic of Korea",
        "Taiwan, Province of China" => "Taiwan",
        "Turkey" => "Türkiye",
    )
    df_mer[!, :location_old] = df_mer[:, :location]
    df_mer[!, :location] = replace(df_mer[:, :location], rep_dic...)
    return df_mer
end

function create_imp_prob_dataframe(df_imp1, df_imp2, df_imp3, df_imp4;
    sort_col::Symbol=:none
)
    df_mer = outerjoin(df_imp1, df_imp2, on=[:iso_code, :location])
    df_mer = outerjoin(df_mer, df_imp3, on=[:iso_code, :location])
    df_mer = outerjoin(df_mer, df_imp4, on=[:iso_code, :location])
    return create_imp_prob_dataframe(df_mer; sort_col=sort_col)
end

function insert_colormap_info!(df_mer)
    # Colormap settings
    cmap = palette(:default)
    blue, green, red = cmap[1], cmap[3], cmap[7]
    rep_dic = Dict(
        "Western Asia" => red,
        "Southern Asia" => red,
        "South-eastern Asia" => blue,
        "Eastern Asia" => green,
        "Central Asia" => red,
    )
    df_mer[!, :color] = replace(df_mer[:, :Region], rep_dic...)
    return df_mer
end

function insert_obs_date(df_mer)
    df_obs = read_observed_imp_data()
    df_mer = leftjoin(
        df_mer,
        df_obs[:, [:iso_code, :first_repo]], on=:iso_code
    )
    df_mer[:, :first_repo] = coalesce.(
        df_mer[:, :first_repo], "2025-12-31"
    ) .|> Date
    return df_mer
end

function insert_tick_labels!(df_mer)
    sort!(df_mer, [:prob1])
    yticks_label = df_mer[:, :location]

    df_obs = read_observed_imp_data()
    df_mer_obs = leftjoin(
        df_mer,
        df_obs[:, [:iso_code, :asia_imp_flag, :obs_imp_flag]],
        on=:iso_code)
    sort!(df_mer_obs, [:prob1])


    cond_double = coalesce.(df_mer_obs[:, :asia_imp_flag], 0) .== 1
    double_ast_cnts = df_mer_obs[cond_double, :location]

    cond_single = coalesce.(df_mer_obs[:, :obs_imp_flag], 0) .== 1
    cond_single = cond_single .& (cond_double .== 0)
    single_ast_cnts = df_mer_obs[cond_single, :location]

    for cnt in single_ast_cnts
        yticks_label = replace.(yticks_label, cnt => cnt * "*")
    end
    for cnt in double_ast_cnts
        yticks_label = replace.(yticks_label, cnt => cnt * "**")
    end
    for i in 1:length(yticks_label)
        l = yticks_label[i]
        if endswith(l, "*") == false
            yticks_label[i] *= "  "
        elseif endswith(l, "**") == false
            yticks_label[i] *= " "
        end
    end
    df_mer[:, :tick_label] = yticks_label
    return nothing
end

function filtered_imp_prob(I_inc, rename_col)
    df_mer = summarise_imp_prob(I_inc)
    cols = [:iso_code, :location, :imp_prob]
    df_mer = df_mer[:, cols]
    # Exclude JPN
    df_mer = @pipe filter(x -> x[:iso_code] != "JPN", df_mer) |>
                   DataFrames.rename(_, :imp_prob => rename_col)
    DataFrames.transform!(df_mer, rename_col => (x -> Float64.(x)) => rename_col)
    return df_mer
end

function global_final_size(fs::Vector{<:Real}; ylabel="")
    q025, q050, q950 = @pipe [
        quantile(fs, q) for q in [0.025, 0.50, 0.975]
    ] .|> round(_, digits=1)
    println("Final outbreak size: $q050 ($q025, $q950)")

    fs = @pipe log10.(fs)

    bins = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0]
    #bins = [100, 300, 1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000, 3_000_000] .|> log10
    xticks = bins
    xticks_label = ["10²", "", "10³", "", "10⁴", "", "10⁵", "", "10⁶", "", "10⁷"]
    v = histogram(fs,
        bins=bins,
        xticks=(xticks, xticks_label),
        label="",
        ylabel=ylabel, xlabel="Final outbreak size in Asia",
        xlabelfontsize=labelfontsize,
        ylabelfontsize=labelfontsize,
        color=palette(:default)[6],
    )
    v[1][1][:y] = v[1][1][:y] ./ length(fs) * 100
    v = plot(v, ylim=[0, 100])
    return v
end

function global_final_size(I_inc::Array{<:Real,3}; ylabel="")
    fs = sum(I_inc, dims=[2, 3])
    return global_final_size(fs[:, 1, 1]; ylabel=ylabel)
end

"""Histogram for the number of countries experiencing
more than `cutoff` cumulative incidence.
"""
function histogram_I_cum_countries(I_inc::Array{Int64,3}; cutoff=10, ylabel="")
    n_sim = size(I_inc)[1]
    n_country = size(I_inc)[3]
    I_cum = cumsum(I_inc, dims=2)
    I_bool = I_cum[:, end, :] .> cutoff
    n_out_cnts = sum(I_bool, dims=2)[:, 1] .- 1
    bins = [0, 1, 5, 10, 20, 30, 40]
    p = histogram(n_out_cnts,
        bins=bins,
        xlabel="Number of countries\nwith importation",
        ylabel=ylabel,
        fmt=:png,
        label="",
        xlabelfontsize=labelfontsize,
        ylabelfontsize=labelfontsize,
        color=palette(:default)[6],
    )
    p[1][1][:y] = p[1][1][:y] ./ (n_sim - 1) * 100
    p = plot(p, ylim=[0, 100])
end

function barplot_I_cum_countries(n_cnts::Vector; ylabel="")
    q025, q050, q950 = @pipe [
        quantile(n_cnts, q) for q in [0.025, 0.50, 0.975]
    ] .|> round(_, digits=1)
    println("Number of countries experiencing importation: $q050 ($q025, $q950)")

    n_sim = length(n_cnts)
    n_out_cnts = n_cnts .- 1
    labels = ["0", "1-5", "6-10", "11-20", ">20"]
    counts = cut(n_out_cnts, [0, 1, 6, 11, 21, 42],
        labels=labels
    ) |> countmap
    y = [get(counts, l, 0) / n_sim * 100 for l in labels]
    pl = bar(labels, y,
        label="",
        ylabel=ylabel,
        xlabel="Number of countries \nwhere  ≥1 importations",
        xrotation=45,
        ylim=[0, 100],
        xlabelfontsize=labelfontsize,
        ylabelfontsize=labelfontsize,
        color=palette(:default)[6],
    )
end

function barplot_I_cum_countries(I_inc::Array{Float64,3};
    cutoff=0, ylabel="")
    n_sim = size(I_inc)[1]
    n_country = size(I_inc)[3]
    I_cum = cumsum(I_inc, dims=2)
    I_bool = I_cum[:, end, :] .> cutoff
    n_cnts = sum(I_bool, dims=2)[:, 1]
    return barplot_I_cum_countries(n_cnts; ylabel=ylabel)
end


"""
    filter_I_inc(I_inc; cut_off=10, tp="Inc_num")

- tp::String: Takes the following params
    - "Inc_num": Conditioned by number of infections outside Japan.
    - "Chnia" (index 8): Conditioned on the number of infection in China.
    - "Korea" (index 18): Conditioned on the number of infection in Korea.
    - "Taiwan" (index 35): Conditioned on the number of infection in Taiwan.
"""
function filter_I_inc(I_inc; cut_off=100, tp="Inc_num")
    cond = filter_I_inc_cond(I_inc; cut_off=cut_off, tp=tp)
    return I_inc[cond, :, :]
end

"""
    filter_I_inc_cond

# Arguments
- `fil_week`: Used for filtering period.
    `fil_week` * 7 is used for maximum index.
    24 is the lenngth of the target value.
"""
function filter_I_inc_cond(I_inc::Array{Float64,3};
    cut_off=100, tp="Inc_num", fil_week=24)
    if tp == "Inc_num"
        fs = @pipe I_inc[:, begin:(fil_week*7), :] |>
                   sum(_, dims=[2, 3])[:, 1, 1]
        jpn_fs = sum(I_inc[:, :, ind0_cnt], dims=2) |> (x -> x[:, 1])
        fs_non_jpn = fs .- jpn_fs
        cond = fs_non_jpn .>= cut_off
    else
        if tp == "China"
            ind = 8
        elseif tp == "Korea"
            ind = 18
        elseif tp == "Taiwan"
            ind = 35
        else
            return nothing
        end
        cnt_fs = sum(I_inc[:, begin:(fil_week*7), ind], dims=2) |> (x -> x[:, 1])
        cond = cnt_fs .>= cut_off
    end
    return cond
end

function check_cond_single_I_inc(I_inc::Matrix{Int64};
    tp="Korea", cut_off=10, fil_week=24
)::Bool
    if tp == "China"
        ind = 8
    elseif tp == "Korea"
        ind = 18
    elseif tp == "Taiwan"
        ind = 35
    else
        return nothing
    end
    cnt_fs = sum(I_inc[begin:(fil_week*7), ind])
    cond = cnt_fs >= cut_off
    return cond
end

"""
...
Note:
    `sim_jpn.jld2` should be moved by yourself.
...
"""
function move_eligible_file_to_directory(vars, path_to; tp="Korea")
    if isdir(path_to) == true
        error("To continue, delete existing $(path_to)")
    else
        mkdir(path_to)
    end
    path_lis = fetch_sim_paths(vars.path_sim_res)
    @showprogress for path in path_lis
        res = JLD2.load(path)["res"]
        I_inc = sum(res.I_new, dims=[3])[:, :, 1]
        flag = check_cond_single_I_inc(I_inc; tp="Korea")
        if flag == true
            file_name = basename(path)
            cp(path, "$(path_to)/$(file_name)")
        end
    end
end

create_record_obj() = return (
    fs=[], imp_flag=[], num_cnt=[],
    jpn_weekly=[],
)
function append_rec_info!(r_obj::NamedTuple, I_inc::Matrix)
    fs = sum(I_inc)
    imp_flag = sum(I_inc, dims=1)[1, :] .>= 1
    num_cnt = sum(imp_flag)
    # To apply `one_country_I_inc_to_weekly`, this should be matrix.
    I_inc = I_inc[:, [15]] |> transpose |> Matrix
    jpn_weekly = one_country_I_inc_to_weekly(I_inc)[1, :]
    append!(r_obj.fs, fs)
    append!(r_obj.imp_flag, [imp_flag])
    append!(r_obj.num_cnt, num_cnt)
    append!(r_obj.jpn_weekly, [jpn_weekly])
end

function record_summary_statistics(vars, path_rec; nmax=100_000)
    r_obj = create_record_obj()
    Korea_ind = []
    Taiwan_ind = []
    China_ind = []
    path_lis = fetch_sim_paths(vars.path_sim_res)
    nmax = minimum([nmax, length(path_lis)])
    @showprogress for (i, path) in enumerate(path_lis[1:nmax])
        res = JLD2.load(path)["res"]
        I_inc = sum(res.I_new, dims=[3])[:, :, 1]
        append_rec_info!(r_obj, I_inc)

        flag_Korea = check_cond_single_I_inc(I_inc; tp="Korea")
        append!(Korea_ind, flag_Korea)

        flag_Taiwan = check_cond_single_I_inc(I_inc; tp="Taiwan")
        append!(Taiwan_ind, flag_Taiwan)

        flag_China = check_cond_single_I_inc(I_inc; tp="China")
        append!(China_ind, flag_China)

    end
    JLD2.jldsave(path_rec,
        r_obj=r_obj, Korea_ind=Korea_ind, Taiwan_ind=Taiwan_ind, China_ind=China_ind)
end

function calc_imp_prob(r_obj::NamedTuple; ind_num=[])
    imp_flag = r_obj.imp_flag
    if length(ind_num) != 0 # Filter
        imp_flag = imp_flag[ind_num]
    end

    n_sim = length(imp_flag)
    imp_prob = fill(0, 42)
    for i in 1:n_sim
        imp_prob .+= imp_flag[i]
    end
    imp_prob = imp_prob ./ n_sim
    return imp_prob .|> Float64
end

function add_hatched_pattern(pl, prob)
    ind = argmax(prob .== 1.0)
    x = [0, prob[ind] * 100]
    y = [ind, ind] .- 0.40
    plot!(pl, x, y, fillrange=y .+ 0.85,
        fillstyle=:/, label=false, color=:black
    )
end

function visualise_imp_prob_global_fs_num_cnts(
    I_inc1, I_inc2, I_inc3, I_inc4,
    title1, title2, title3, title4;
    sort_col=:prob1, kwds=nothing
)
    kwds = deepcopy(kwds)
    println(size.([I_inc1, I_inc2, I_inc3, I_inc4]))
    # Prepare merged dataframes
    df_imp1 = filtered_imp_prob(I_inc1, :prob1)
    df_imp2 = filtered_imp_prob(I_inc2, :prob2)
    df_imp3 = filtered_imp_prob(I_inc3, :prob3)
    df_imp4 = filtered_imp_prob(I_inc4, :prob4)

    df_mer = create_imp_prob_dataframe(df_imp1, df_imp2, df_imp3, df_imp4;
        sort_col=sort_col)
    insert_colormap_info!(df_mer)
    df_mer = insert_obs_date(df_mer)
    sort!(df_mer, [sort_col], rev=false)
    insert_tick_labels!(df_mer)
    kwds[:color] = df_mer[:, :color]
    colors = df_mer[:, :color]

    pl1 = bar(
        df_mer[:, :prob1] * 100,
        yticks=(1:n_cnt, df_mer[:, :tick_label]),
        #title="natsal 1 year\nip 7 days";
        title=title1;
        kwds...
    )
    add_hatched_pattern(pl1, df_mer[:, :prob1])

    pl2 = bar(
        df_mer[:, :prob2] * 100,
        yticks=yticks_empty,
        title=title2;
        kwds...
    )
    add_hatched_pattern(pl2, df_mer[:, :prob2])

    pl3 = bar(
        df_mer[:, :prob3] * 100,
        yticks=yticks_empty,
        title=title3;
        kwds...
    )
    add_hatched_pattern(pl3, df_mer[:, :prob3])

    pl4 = bar(
        df_mer[:, :prob4] * 100,
        yticks=yticks_empty,
        title=title4;
        kwds...
    )
    add_hatched_pattern(pl4, df_mer[:, :prob4])

    # Legened
    cmap = palette(:default)
    blue, green, red = cmap[1], cmap[3], cmap[7]
    bar!(pl1, [0 0 0], #[1 1 1],
        legend=(0.3, 0.2),
        labels=[" Eastern Asia" " South-eastern Asia" " Other Asia"],
        color=[blue green red],
        lw=1,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        orientation=:horizontal,
    )

    # global final size epidemics,
    pl1_fs = global_final_size(I_inc1; ylabel="Proportion (%)")
    pl2_fs = global_final_size(I_inc2)
    pl3_fs = global_final_size(I_inc3)
    pl4_fs = global_final_size(I_inc4)

    # Imp countries
    pl1_imp = barplot_I_cum_countries(I_inc1, cutoff=0, ylabel="Proportion (%)")
    pl2_imp = barplot_I_cum_countries(I_inc2, cutoff=0)

    pl3_imp = barplot_I_cum_countries(I_inc3, cutoff=0)
    pl4_imp = barplot_I_cum_countries(I_inc4, cutoff=0)
    annotate!(pl1, (-0.7, 1.03), text("A", :black, :left, 20))
    annotate!(pl1_fs, (-0.7, 0.9), text("B", :black, :left, 20))
    annotate!(pl1_imp, (-0.7, 0.9), text("C", :black, :left, 20))

    # adjust figures
    pl = [
        pl1, pl2, pl3, pl4,
        pl1_fs, pl2_fs, pl3_fs, pl4_fs,
        pl1_imp, pl2_imp, pl3_imp, pl4_imp
    ]
    l = @layout [
        a{0.8h} b c d;
        e f g h;
        i j k l
    ]
    plot(pl..., layout=l, size=(900, 1000), fmt=:png, dpi=300)
end

"""
...
Args
- `path_rec`: Recording path for summarised statistics.
...
"""
function visualise_imp_prob_global_fs_num_cnts_large_data(
    path_rec, title1, title2, title3, title4;
    sort_col=:prob1, kwds=nothing
)
    sum_obj = load(path_rec)
    r_obj = sum_obj["r_obj"]

    Korea_ind = [i for (i, x) in enumerate(sum_obj["Korea_ind"])
                 if x == true]
    Taiwan_ind = [i for (i, x) in enumerate(sum_obj["Taiwan_ind"])
                  if x == true]
    China_ind = [i for (i, x) in enumerate(sum_obj["China_ind"])
                 if x == true]

    # Importation prob.
    df_mer = CSV.read(path_pop, DataFrame)
    df_mer[!, :prob1] = calc_imp_prob(r_obj)
    df_mer[!, :prob2] = calc_imp_prob(r_obj; ind_num=Korea_ind)
    df_mer[!, :prob3] = calc_imp_prob(r_obj; ind_num=Taiwan_ind)
    df_mer[!, :prob4] = calc_imp_prob(r_obj; ind_num=China_ind)
    filter!(x -> x["iso_code"] ≠ "JPN", df_mer)
    df_mer = create_imp_prob_dataframe(df_mer; sort_col=sort_col)

    # Edit dataframe for visualisation.
    insert_colormap_info!(df_mer)
    df_mer = insert_obs_date(df_mer)
    sort!(df_mer, [sort_col], rev=false)
    insert_tick_labels!(df_mer)
    kwds[:color] = df_mer[:, :color]
    colors = df_mer[:, :color]

    n_cnt = nrow(df_mer)
    pl1 = bar(
        df_mer[:, :prob1] * 100,
        yticks=(1:n_cnt, df_mer[:, :tick_label]),
        #title="natsal 1 year\nip 7 days";
        title=title1;
        kwds...
    )
    add_hatched_pattern(pl1, df_mer[:, :prob1])

    pl2 = bar(
        df_mer[:, :prob2] * 100,
        yticks=yticks_empty,
        title=title2;
        kwds...
    )
    add_hatched_pattern(pl2, df_mer[:, :prob2])

    pl3 = bar(
        df_mer[:, :prob3] * 100,
        yticks=yticks_empty,
        title=title3;
        kwds...
    )
    add_hatched_pattern(pl3, df_mer[:, :prob3])

    pl4 = bar(
        df_mer[:, :prob4] * 100,
        yticks=yticks_empty,
        title=title4;
        kwds...
    )
    add_hatched_pattern(pl4, df_mer[:, :prob4])

    # Legened
    cmap = palette(:default)
    blue, green, red = cmap[1], cmap[3], cmap[7]
    bar!(pl1, [0 0 0], #[1 1 1],
        legend=(0.3, 0.2),
        labels=[" Eastern Asia" " South-eastern Asia" " Other Asia"],
        color=[blue green red],
        # color=[palette(:algae)[100] palette(:ice)[190] palette(:amp)[150]],
        lw=1,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        orientation=:horizontal,
    )

    # global final size epidemics,
    fs = r_obj.fs .|> Float64
    pl1_fs = global_final_size(fs; ylabel="Proportion (%)")
    pl2_fs = global_final_size(fs[Korea_ind])

    pl3_fs = global_final_size(fs[Taiwan_ind])
    pl4_fs = global_final_size(fs[China_ind])

    # Imp countries
    pl1_imp = barplot_I_cum_countries(r_obj.num_cnt, ylabel="Proportion (%)")
    pl2_imp = barplot_I_cum_countries(r_obj.num_cnt[Korea_ind])

    pl3_imp = barplot_I_cum_countries(r_obj.num_cnt[Taiwan_ind])
    pl4_imp = barplot_I_cum_countries(r_obj.num_cnt[China_ind])

    annotate!(pl1, (-0.7, 1.03), text("A", :black, :left, 20))
    annotate!(pl1_fs, (-0.7, 0.9), text("B", :black, :left, 20))
    annotate!(pl1_imp, (-0.7, 0.9), text("C", :black, :left, 20))

    # adjust figures
    pl = [
        pl1, pl2, pl3, pl4,
        pl1_fs, pl2_fs, pl3_fs, pl4_fs,
        pl1_imp, pl2_imp, pl3_imp, pl4_imp
    ]
    l = @layout [
        a{0.8h} b c d;
        e f g h;
        i j k l
    ]
    plot(pl..., layout=l, size=(900, 1000), fmt=:png, dpi=300)
end

function one_country_I_inc_to_weekly(I_inc::Matrix)
    n_row, n_col = size(I_inc)
    n_point = Int64(n_col / 7 |> floor)
    I_inc_weekly = fill(0, n_row, n_point)
    for i in 1:n_row
        for j in 1:n_point
            st = 1 + (j - 1) * 7
            fin = j * 7
            I_inc_weekly[i, j] = sum(I_inc[i, st:fin])
        end
    end
    return I_inc_weekly
end

function imp_date_vis(imp_dates::DataFrame, df_mer::DataFrame; yticks=true, title="",
    label=true,
)
    cols = [:q05, :q50, :q95]
    df_imp = copy(imp_dates)
    df_imp[:, cols] = ifelse.(isnan.(df_imp[:, cols]), 1000, df_imp[:, cols])
    # For reporting delay including 3 days incubation period and
    # 9 days onset to medical attendance.
    df_imp[:, cols] .+= 12
    df_imp = leftjoin(df_mer, df_imp, on=:location_old => :country)
    insert_tick_labels!(df_imp)
    sort!(df_imp, [:prob1])

    n_country = size(df_imp)[1]
    st_date = Date(2023, 1, 16)
    x_min = Date(2022, 12, 25)
    x_max = Date(2024, 12, 31)

    q50s = (df_imp[:, :q50])
    q50s_date = st_date .+ Day.(q50s)
    xticks = [st_date + Month(i - 1) for i in 1:24]
    yticks = yticks == true ? (1:n_country, df_imp[:, :tick_label]) : (1:n_country, "")


    pl = plot(
        xlim=[x_min.instant.periods.value,
            x_max.instant.periods.value],
        ylim=[0, n_country + 1],
        xrotation=90,
        xticks=(xticks, xticks),
        yticks=yticks,
        #size=(400, 800),
        xlabel="First importation date since 2023",
        title=title, #  * ", Median (5th, 95th)",
        #legend=(0.2,0.2),
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        size=(700, 900),
    )
    label1 = "With ≥1 importations"
    label1_non = "Without importations"
    label2 = "Observed"
    label3 = "No observation by 2023"
    if label == false
        label1 = :none
        label1_non = :none
        label2 = :none
        label3 = :none
    end
    blue, orange, green = palette(:default)[[1, 2, 3]]

    # Draw credible intervals.
    ys = 1:n_country
    dy = 0.25
    x_l = @pipe st_date .+ Day.(df_imp[:, :q05]) |> reshape(_, 1, :)
    x_h = @pipe st_date .+ Day.(df_imp[:, :q95]) |> reshape(_, 1, :)


    label1_lis = vcat(["", "", "", "", "",], label1, label1_non, ["" for i in 8:n_country])
    cond = df_imp[:, :first_repo] .== Date(2025, 12, 31)
    color = @pipe [f == true ? blue : :blue for f in cond] |> reshape(_, 1, :)
    plot!(pl, [x_l; x_h], [ys'; ys'],
        color=color, linealpha=0.3, lw=8,
        label=reshape(label1_lis, 1, :),
    )

    # Draw median marker
    x_del = Day(0)
    x_med = @pipe q50s_date |> reshape(_, 1, :)
    x_l_med = @pipe q50s_date .- x_del |> reshape(_, 1, :)
    x_h_med = @pipe q50s_date .+ x_del |> reshape(_, 1, :)
    dy = 0.1
    plot!(pl, [x_med; x_med], [(ys .- dy)'; (ys .+ dy)'],
        color=color, linealpha=0.4, lw=3.5,
        #markershape=:vline, markersize=4, markeralpha=0.3,
        label="",
    )

    # Draw observed dates.
    cond = df_imp[:, :first_repo] .== Date(2025, 12, 31)
    n_non = sum(cond)
    n_obs = n_country - n_non
    x_l_obs = @pipe [x_min for _ in 1:n_obs] |> reshape(_, 1, :)
    x_h_obs = @pipe df_imp[cond.==false, :first_repo] |> reshape(_, 1, :)
    ys_obs = ys[cond.==false]
    scatter!(pl, df_imp[:, :first_repo], 1:n_country,
        label=label2, color=orange, marker=:circle,
        markerstrokewidth=0.6, markersize=4,
    )

    # Draw non-observed data for cross.
    x_l_cond = @pipe [x_min for _ in 1:n_non] |> reshape(_, 1, :)
    x_h_cond = @pipe [Date(2024, 1, 1) for _ in 1:n_non] |> reshape(_, 1, :)
    ys_cond = ys[cond]
    vspan!(pl, [Date(2024, 1, 1), Date(2025, 12, 31)], color=:black, alpha=0.1, label=:none)
    return pl
end

function save_inc_imp(vars::VisVars)
    I_inc = get_metrix_from_each_sim(vars.path_sim_res, get_incidence_I; nmax=nmax)
    JLD2.save_object(vars.path_inc, I_inc)
    df_imp_prob = summarise_imp_prob(I_inc)
    CSV.write(vars.path_imp, df_imp_prob)
end

function save_filtered_imp(vars::VisVars)
    load_I_inc!(vars; check=true)
    I_inc = filter_I_inc(vars.I_inc; tp="Korea", cut_off=10)
    df_imp_prob = summarise_imp_prob(I_inc)
    CSV.write(vars.path_imp_fil, df_imp_prob)
end

function gen_based_Sankey_diagram(path; nmax=10000)
    path_objs = fetch_sim_paths(path)
    df_sum = DataFrame()

    n = minimum([length(path_objs), nmax])
    @showprogress for i in 1:n
        res = JLD2.load(path_objs[i])["res"]
        global df = DataFrame(res.import_event)
        if size(df)[1] == 0
            continue
        end

        df[!, :sim_index] .= i
        df[!, :export] = replace(df[:, :export_cntry], country_dict...)
        df[!, :import] = replace(df[:, :import_cntry], country_dict...)
        df[!, :ex_im] = df[!, :export] .* ": " .* df[!, :import]
        cond = .!nonunique(df[:, [:sim_index, :export, :import]])
        df = df[cond, :]

        df[!, :gen_index] .= 100000
        df[!, :seq_index] .= 100000

        sort!(df, [:time])
        n_row = size(df)[1]
        seq_index = 1
        pre_df = DataFrame()
        for i in 1:n_row
            if df[i, :export] == "Japan"
                df[i, :gen_index] = 1
                df[i, :seq_index] = seq_index
                seq_index += 1
                pre_df = vcat(pre_df, df[i, :] |> DataFrame)
            else
                cond = df[i, :import] .== pre_df[:, :import]
                if any(cond)
                    continue
                end

                r = filter(x -> x["import"] == df[i, :export], pre_df)[1, :]
                df[i, :gen_index] = r.gen_index + 1
                df[i, :seq_index] = r.seq_index
                pre_df = vcat(pre_df, df[i, :] |> DataFrame)
            end
        end
        df_sum = vcat(df_sum, pre_df)
    end
    return df_sum
end

# Filter datasets
function filter_exp_gen(vars; cut_off=10, tp="Korea")
    df_upd = CSV.read(vars.path_exp_imp_gen, DataFrame)
    load_I_inc!(vars)
    cond = filter_I_inc_cond(vars.I_inc; cut_off=cut_off, tp=tp)
    println("Number of sims: ", sum(cond))
    ind_lis = [i for i in 1:length(cond)][cond]
    ind_fil = intersect(ind_lis, df_upd[:, :sim_index])
    df_fil = filter(x -> x[:sim_index] in ind_fil, df_upd)
    return df_fil
end

"""
...
Args
- `n_sim`: Actual number of simulations.
    After conditioned, the file number is not matched with simulation number.
...
"""
function conditional_beta_SAR(vars::VisVars, cut_off, tp; n_sim=5000)
    load_I_inc!(vars; check=true)
    I_inc = vars.I_inc
    cond = filter_I_inc_cond(I_inc; cut_off=cut_off, tp=tp)

    # Prepare the simulated beta (whether or not values were saved or not.)
    path_objs = fetch_sim_paths(vars.path_sim_res)
    βs = [JLD2.load(path, "β") for path in path_objs]
    println("Number of valid samples: ", length(βs))

    βs_qs = quantile(βs, [0.025, 0.5, 0.975])
    βs_r = round.(βs_qs, digits=2)
    println("β : $(βs_r[2]) ($(βs_r[1]), $(βs_r[3]))")

    SAR_qs = @pipe (1 .- exp.(-βs_qs)) .* 100 .|> round(_, digits=3)
    println("SAR : $(SAR_qs[2]) ($(SAR_qs[1]), $(SAR_qs[3]))")
end

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
    vspan!(pl, [-10, 23], color=:black, alpha=0.2, label=:none)
    for p in ps
        lw = p == 0.5 ? 2 : 1
        ls = p == 0.95 ? :dot : :solid

        q = q_dict[p]
        p_lab = Int64(p * 100)
        plot!(pl, 1:n_week, q .+ 1,
            lw=lw, ls=ls, color="blue",
            label="$(p_lab)th",
            legend_title="Quantile",
        )
    end
    pl
end