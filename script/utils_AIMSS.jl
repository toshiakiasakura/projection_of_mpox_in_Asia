"""Clean dataset in the way written in the paper.
- Sin How Lim 2012 AIDS Behav
https://pmc.ncbi.nlm.nih.gov/articles/PMC4405155/pdf/nihms448989.pdf
"""
function clean_process1(df::DataFrame)::DataFrame
    cond_gender = in.(df[:, :gender] , Ref(["M", "MTF"]))
    cond_age = df[:, :year_birth] .<= 2010 - 18
    cond_nfreq = df[:, :nmsex] .!= "NONE"
    df_fil = df[cond_age .& cond_nfreq .& cond_gender, :]
    println(size(df_fil), ": include `M`, aged >= 18, nfreq >= 1")
    df_fil_tmp = @subset df_fil begin
        :sex_regular .== "N"
        :sex_casual .== "N"
        :sex_comm .== "N"
    end
    df_fil2 = filter(row ->  !(row.responseid in df_fil_tmp.responseid), df_fil)
    println(size(df_fil2), ": exclude participants answering No regular, casual, comm partners.")

    return df_fil2
end

"""Convert categorical values of sexual partner numbers into numerical values.
"""
function rep_and_parse_values(df)
    rep_dic = Dict(" " => 0,  "ONE" => 1,  "2TO5" => 2, "6TO10" => 3,  "11TO50" => 4,  "50PLUS" => 5)
    # Remove contradicted answers.
    df_tmp = @chain df_fil2 begin
        @transform(
            :nmsex_num = replace.(:nmsex, rep_dic...),
            :nregular_num = replace.(:nregular, rep_dic...),
            :ncasual_num = replace.(:ncasual, rep_dic...),
            :ncomm_num = replace.(:ncomm, rep_dic...)
        )
    end
    lis_conv = [:nmsex_num, :nregular_num, :ncasual_num, :ncomm_num]
    df_tmp[!, lis_conv] = parse.(Int64, df_tmp[:, lis_conv])
    return df_tmp
end

function crosstab(df::DataFrame, cols)
    tab = combine(groupby(df, cols), nrow => :Count)
    tab = @pipe unstack(tab, :nmsex_num, :Count) .|> coalesce(_, 0)
    # Add "All" row
    v = sum.(eachcol(tab[:, 2:end]))
    push!(tab, ["All", v...])
    # Add "All" column
    tab[!, "All"] = sum(eachcol(tab[:, 2:end]))
    return tab
end

function convert2longformat(tab::DataFrame)::DataFrame
    n_col = ncol(tab)
    tab_long = stack(tab[:, 1:(n_col-1)], 2:(n_col-1))
    tab_long = leftjoin(tab_long, tab[:, [:countryname, :All]], on=:countryname)
    order = [ "All", "China (mainland)", "Singapore", "Malaysia", "Taiwan", "Hong Kong",
            "Thailand", "Japan", "Indonesia", "Philippines", "Vietnam", "Republic of Korea" ]
    v = CategoricalArray(tab_long[:, :countryname])
    levels!(v, order)
    tab_long[:, :countryname_order] = v
    return tab_long
end

"""yerr is given by Quesenberry-Hurst method.
"""
function add_yerr!(tab_long::DataFrame)::Nothing
    tab_long[!, :yerr_l] .= NaN
    tab_long[!, :yerr_h] .= NaN
    sort!(tab_long, :variable)
    for g in groupby(tab_long, :countryname)
        v = g[:, :value]
        prop = g[:, :prop]
        conf = confint(PowerDivergenceTest(v); level=0.95, method=:quesenberry_hurst)
        g[:, :yerr_l] = prop .- [c[1] .* 100 for c in conf]
        g[:, :yerr_h] = [c[2] .* 100 for c in conf] .- prop
    end
end

function get_pdf_4w_1y()
    lower = 0.5
    bins = [lower, 2, 6, 11, 51, Inf]
    period_y = 0.5
    α, κ = NATSAL4W_PARMS
    wb_4w = truncated(WeibullPareto(α, κ * (1 / period_y)^α); lower = lower)
    α, κ = NATSAL1Y_PARMS
    wb_1y = truncated(WeibullPareto(α, κ * (1 / period_y)^α); lower = lower)

    p_pdf_4w = @pipe cdf.(wb_4w, bins) |> (x -> x[2:end] .- x[1:(end-1)])
    p_pdf_1y = @pipe cdf.(wb_1y, bins) |> (x -> x[2:end] .- x[1:(end-1)])
    normalize!(p_pdf_4w, 1)
    normalize!(p_pdf_1y, 1)
    return (p_pdf_4w, p_pdf_1y)
end

"""
Note: rows in tab are given in an order which is used for the visualisation.
"""
function add_sample_size(tab, tab_long)
    tab = @chain tab begin
        @transform :cname = :countryname .* " (n=" .* string.( :All) .* ")"
    end
    v = tab[:, :cname]
    v_cate = CategoricalArray(v)
    levels!(v_cate, v)
    tab[:, :cname_cate] = v_cate
    tab_long = leftjoin(tab_long, tab[:, [:countryname, :cname_cate]], on=:countryname)
    return tab_long
end

function add_natsal_data!(tab_long::DataFrame; exclude_one=false)::Nothing
    p_pdf_4w, p_pdf_1y = get_pdf_4w_1y()

    if exclude_one == true
        p_pdf_4w = normalize(p_pdf_4w[2:5], 1)
        p_pdf_1y = normalize(p_pdf_1y[2:5], 1)
    end
    adj = exclude_one == true ? 1 : 0

    for (i, p) in enumerate(p_pdf_4w)
        cname = "Natsal 4-week"
        push!(tab_long, [cname, string(i + adj), 0, 0, cname, p*100, NaN, NaN, cname])
    end
    for (i, p) in enumerate(p_pdf_1y)
        cname = "Natsal 1-year"
        push!(tab_long, [cname, string(i + adj), 0, 0, cname, p*100, NaN, NaN, cname])
    end
end

function prepare_tab_long1(df_ana)
    cols = [:countryname, :nmsex_num]
    tab = crosstab(df_ana, cols)
    sort!(tab, :All, rev=false)
    tab = @subset tab :All.> 100
    tab_long = convert2longformat(tab)
    tab_long = @transform(tab_long, :prop = :value./:All.* 100)
    add_yerr!(tab_long)
    tab_long = add_sample_size(tab, tab_long)
    add_natsal_data!(tab_long)
    return tab_long
end

function prepare_tab_long2(df_ana)
    cols = [:countryname, :nmsex_num]
    tab_ori = crosstab(df_ana, cols)
    tab_ori = @rename(tab_ori, :All_ori = :All)

    df_ana_fil = @subset(df_ana, :nmsex_num .!= 1)
    tab = crosstab(df_ana_fil, cols)
    tab = leftjoin(tab, tab_ori[:, [:countryname, :All_ori]], on=:countryname)

    # Sort based on the total number to match tab_long1.
    sort!(tab, :All_ori, rev=false)
    tab = tab[:, begin:(end-1)]
    tab = @subset tab :All.> 80
    tab_long = convert2longformat(tab)
    tab_long = @transform(tab_long, :prop = :value./:All.* 100)
    add_yerr!(tab_long)
    tab_long = add_sample_size(tab, tab_long)
    add_natsal_data!(tab_long; exclude_one=true)
    return tab_long
end

function plot_AIMSS_bar1(tab_long, kwds_bar)
    xticks = ["1", "2-5", "6-10", "11-50", ">50"]
    t = tab_long
    pl = groupedbar(t[:, :variable], t[:, :prop]; group=t[:, :cname_cate],
            yerr=(t[:, :yerr_l], t[:, :yerr_h]), xticks=(0.5:1:5, xticks),
            legend=(0.75, 0.95),
            kwds_bar...
    )
    pl
end

function plot_AIMSS_bar2(tab_long, kwds_bar)
    xticks = ["2-5", "6-10", "11-50", ">50"]
    t = tab_long
    pl = groupedbar(t[:, :variable], t[:, :prop]; group=t[:, :cname_cate],
            yerr=(t[:, :yerr_l], t[:, :yerr_h]), xticks=(0.5:1:5, xticks),
            legend= true,
            kwds_bar...
    )
    return pl
end

function steplines!(pl::Plots.Plot, v::Vector, label; color = palette(:default)[1])
    for (i,p) in enumerate(v)
        label = i == 1 ? label : nothing
        plot!(pl, [i-0.97,i-0.03], [p, p].*100, color=color, label=label, lw=1.2)
    end
end