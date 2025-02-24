# Visualisation of intercountry simulation results.

include("data_path.jl")
include("./model_intercountry.jl")

N_CNT_ECL_JPN = 41
BAR_LABEL_FONT_SIZE = 8
BAR_KWDS = Dict(
	:xlim => [0, 105],
	:orientation => :horizontal,
	:ylim => [0, N_CNT_ECL_JPN + 1],
	:label => "",
	:left_margin => 7Plots.pt,
	:right_margin => 0Plots.pt,
	:upper_margin => 40Plots.pt,
	:xlabel => "Simulated importation\nprobability (%)",
	:titlefontsize => 10,
	:ytickfontsize => BAR_LABEL_FONT_SIZE,
	:xlabelfontsize => BAR_LABEL_FONT_SIZE,
	:ylabelfontsize => BAR_LABEL_FONT_SIZE,
)
CNT_REP_DIC = Dict(
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

function get_metrix_from_each_sim(path::String, f; nmax::Int64 = 100_000)
	path_objs = fetch_sim_paths(path)
	n = length(path_objs)
	n = minimum([n, nmax])
	sim_res::Array{Any, 1} = fill(Inf, n)

	ThreadsX.foreach(1:n) do (i)
		try
			res = JLD2.load(path_objs[i])["res"]
			sim_res[i] = f(res)
		catch
			error("Index $i")
		end
	end
	sim_mat = vec_matrix_to_matrix(sim_res)
	return sim_mat
end

function fetch_initial_imports(res::ResultInterCountrySEIRModel)
	df_imp = DataFrame(res.import_event)
	initial_imports = fill(0, res.n_country, res.days)
	if nrow(df_imp) == 0
		return initial_imports
	end

	df_min = @pipe groupby(df_imp, :import_cntry) |> combine(_, :time => minimum => :time)
	for i in 1:res.n_country
		v = filter(x -> x[:import_cntry] == i, df_min)
		if nrow(v) == 1
			time = v[1, :time]
			initial_imports[i, time] += 1
		end
	end
	return initial_imports
end

function load_I_inc!(vars::VisVars; check = false)
	if (check == true) & (isinf(vars.I_inc[1, 1, 1]) == false)
		return nothing
	end
	vars.I_inc = JLD2.load_object(vars.path_inc)
	return nothing
end

get_incidence_I(model::ResultInterCountrySEIRModel) = model.I_new

function create_imp_prob_dataframe(df_mer::DataFrame; sort_col::Symbol = :none)
	# Add Region information.
	df_UN = @select CSV.read("../data/UN_Asia_list.csv", DataFrame) :Code :Region
	df_mer = leftjoin(df_mer, df_UN, on = :iso_code => :Code)

	println("Regions of the following countires were imputed with Eastern Asia")
	@with df_mer begin
		cond = :Region .|> ismissing
		df_mer[cond, :] |> display
		:Region[cond] .= "Eastern Asia"
		:Region |> unique |> display
	end

	df_mer = sort(df_mer, [sort_col], rev = false)
	@transform!(df_mer, :location_old = :location)
	df_mer[!, :location] = replace(df_mer[:, :location], CNT_REP_DIC...)
	return df_mer
end

function create_imp_prob_dataframe(df_imp1, df_imp2, df_imp3, df_imp4;
	sort_col::Symbol = :none,
)
	df_mer = outerjoin(df_imp1, df_imp2, df_imp3, df_imp4, on = [:iso_code, :location])
	return create_imp_prob_dataframe(df_mer; sort_col = sort_col)
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
end

function insert_obs_date(df_mer)
	df_obs = @select read_observed_imp_data() :iso_code :first_repo
	df_mer = leftjoin(df_mer, df_obs, on = :iso_code)
	df_mer[:, :first_repo] = coalesce.(
		df_mer[:, :first_repo], "2025-12-31",
	) .|> Date
	return df_mer
end

function insert_tick_labels!(df_mer)
	yticks_label = df_mer[:, :location]

	df_obs = @select read_observed_imp_data() :iso_code :asia_imp_flag :obs_imp_flag
	df_mer_obs = leftjoin(df_mer, df_obs, on = :iso_code)

	cond_double = coalesce.(df_mer_obs[:, :asia_imp_flag], 0) .== 1
	ast2_cnts = df_mer_obs[cond_double, :location]
	ast1_cnts = @subset(df_mer_obs,
		(:obs_imp_flag .== 1) .& (cond_double .== 0)
	)[:, :location]

	cond2 = in.(yticks_label, Ref(ast2_cnts))
	cond1 = in.(yticks_label, Ref(ast1_cnts))
	df_mer[!, :tick_label] = @pipe yticks_label |>
								   ifelse.(cond2, _ .* "**", _) |>
								   ifelse.(cond1, _ .* "* ", _) |>
								   ifelse.(.~(cond1 .| cond2), _ .* "  ", _)
	return nothing
end

function filtered_imp_prob(I_inc::Array{Float64, 3}, rename_col)
	df_mer = summarise_imp_prob(I_inc)
	df_mer = @chain df_mer begin
		@select :iso_code :location :imp_prob
		@subset :iso_code .!= "JPN"
		@rename $rename_col = :imp_prob
		@transform $rename_col = Float64.($rename_col)
	end
	return df_mer
end

function global_final_size(fs::Vector{<:Real}; ylabel = "")
	q025, q050, q950 = @pipe [
		quantile(fs, q) for q in [0.025, 0.50, 0.975]
	] .|> round(_, digits = 1)
	println("Final outbreak size: $q050 ($q025, $q950)")

	fs = @pipe log10.(fs)

	bins = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0]
	xticks = bins
	xticks_label = ["10²", "", "10³", "", "10⁴", "", "10⁵", "", "10⁶", "", "10⁷"]
	v = histogram(fs,
		bins = bins,
		xticks = (xticks, xticks_label),
		label = "",
		ylabel = ylabel, xlabel = "Final outbreak size in Asia",
		xlabelfontsize = BAR_LABEL_FONT_SIZE,
		ylabelfontsize = BAR_LABEL_FONT_SIZE,
		color = palette(:default)[6],
	)
	v[1][1][:y] = v[1][1][:y] ./ length(fs) * 100
	v = plot(v, ylim = [0, 100])
	return v
end

function global_final_size(I_inc::Array{<:Real, 3}; ylabel = "")
	fs = sum(I_inc, dims = [2, 3])
	return global_final_size(fs[:, 1, 1]; ylabel = ylabel)
end

"""Histogram for the number of countries experiencing
more than `cutoff` cumulative incidence.
"""
function histogram_I_cum_countries(I_inc::Array{Int64, 3}; cutoff = 10, ylabel = "")
	n_sim = size(I_inc)[1]
	n_country = size(I_inc)[3]
	I_cum = cumsum(I_inc, dims = 2)
	I_bool = I_cum[:, end, :] .> cutoff
	n_out_cnts = sum(I_bool, dims = 2)[:, 1] .- 1
	bins = [0, 1, 5, 10, 20, 30, 40]
	p = histogram(n_out_cnts,
		bins = bins,
		xlabel = "Number of countries\nwith importation",
		ylabel = ylabel,
		fmt = :png,
		label = "",
		xlabelfontsize = BAR_LABEL_FONT_SIZE,
		ylabelfontsize = BAR_LABEL_FONT_SIZE,
		color = palette(:default)[6],
	)
	p[1][1][:y] = p[1][1][:y] ./ (n_sim - 1) * 100
	p = plot(p, ylim = [0, 100])
end

function barplot_I_cum_countries(n_cnts::Vector; ylabel = "",
	xlabel = "Number of countries \nwhere  ≥1 importations",
	xrotation=45,
)
	q025, q050, q950 = @pipe [
		quantile(n_cnts, q) for q in [0.025, 0.50, 0.975]
	] .|> round(_, digits = 1)
	println("Number of countries experiencing importation: $q050 ($q025, $q950)")

	n_sim = length(n_cnts)
	n_out_cnts = n_cnts .- 1
	labels = ["0", "1-5", "6-10", "11-20", ">20"]
	counts = cut(n_out_cnts, [0, 1, 6, 11, 21, 42],
		labels = labels,
	) |> countmap
	y = [get(counts, l, 0) / n_sim * 100 for l in labels]
	pl = bar(labels, y,
		label = "",
		ylabel = ylabel,
		xlabel = xlabel,
		xrotation = xrotation,
		ylim = [0, 100],
		xlabelfontsize = BAR_LABEL_FONT_SIZE,
		ylabelfontsize = BAR_LABEL_FONT_SIZE,
		color = palette(:default)[6],
	)
end

function barplot_I_cum_countries(I_inc::Array{Float64, 3};
	cutoff = 0, kwds...)
	n_sim = size(I_inc)[1]
	n_country = size(I_inc)[3]
	I_cum = cumsum(I_inc, dims = 2)
	I_bool = I_cum[:, end, :] .> cutoff
	n_cnts = sum(I_bool, dims = 2)[:, 1]
	return barplot_I_cum_countries(n_cnts; kwds...)
end


"""
	filter_I_inc(I_inc; cut_off=10, tp="Inc_num")

- tp::String: Takes the following params
	- "Inc_num": Conditioned by number of infections outside Japan.
	- "Chnia" (index 8): Conditioned on the number of infection in China.
	- "Korea" (index 18): Conditioned on the number of infection in Korea.
	- "Taiwan" (index 35): Conditioned on the number of infection in Taiwan.
"""
function filter_I_inc(I_inc; cut_off = 100, tp = "Inc_num")
	cond = filter_I_inc_cond(I_inc; cut_off = cut_off, tp = tp)
	return I_inc[cond, :, :]
end

"""
	filter_I_inc_cond

# Arguments
- `fil_week`: Used for filtering period.
	`fil_week` * 7 is used for maximum index.
	24 is the lenngth of the target value.

"""
function filter_I_inc_cond(I_inc::Array{Float64, 3};
	cut_off = 100, tp = "Inc_num", fil_week = 24)
	if tp == "Inc_num"
		fs = @pipe I_inc[:, begin:(fil_week*7), :] |>
				   sum(_, dims = [2, 3])[:, 1, 1]
		jpn_fs = sum(I_inc[:, :, ind0_cnt], dims = 2) |> (x -> x[:, 1])
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
		cnt_fs = sum(I_inc[:, begin:(fil_week*7), ind], dims = 2) |> (x -> x[:, 1])
		cond = cnt_fs .>= cut_off
	end
	return cond
end

function check_cond_single_I_inc(I_inc::Matrix{Int64};
	tp = "Korea", cut_off = 10, fil_week = 24,
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
function move_eligible_file_to_directory(vars::VisVars, vars_to::VisVars; tp = "Korea")
	path_to = vars_to.path_sim_res
	if isdir(path_to) == true
		error("To continue, delete existing $(path_to)")
	else
		mkdir(path_to)
	end
	path_lis = fetch_sim_paths(vars.path_sim_res)

	Threads.@threads for path in path_lis
		res = JLD2.load(path)["res"]
		I_inc = res.I_new
		flag = check_cond_single_I_inc(I_inc; tp = "Korea")
		if flag == true
			file_name = basename(path)
			cp(path, "$(path_to)/$(file_name)")
		end
	end
end

create_record_obj() = return (
	fs = [], imp_flag = [], num_cnt = [],
	jpn_weekly = [],
)
function append_rec_info!(r_obj::NamedTuple, I_inc::Matrix)
	fs = sum(I_inc)
	imp_flag = sum(I_inc, dims = 1)[1, :] .>= 1
	num_cnt = sum(imp_flag)
	# To apply `one_country_I_inc_to_weekly`, this should be matrix.
	I_inc = I_inc[:, [15]] |> transpose |> Matrix
	jpn_weekly = one_country_I_inc_to_weekly(I_inc)[1, :]
	append!(r_obj.fs, fs)
	append!(r_obj.imp_flag, [imp_flag])
	append!(r_obj.num_cnt, num_cnt)
	append!(r_obj.jpn_weekly, [jpn_weekly])
end

function record_summary_statistics(vars::VisVars; nmax = 100_000)
	r_obj = create_record_obj()
	Korea_ind = []
	Taiwan_ind = []
	China_ind = []
	path_lis = fetch_sim_paths(vars.path_sim_res)
	nmax = minimum([nmax, length(path_lis)])
	@showprogress for (i, path) in enumerate(path_lis[1:nmax])
		res = JLD2.load(path)["res"]
		I_inc = res.I_new
		append_rec_info!(r_obj, I_inc)

		flag_Korea = check_cond_single_I_inc(I_inc; tp = "Korea")
		append!(Korea_ind, flag_Korea)

		flag_Taiwan = check_cond_single_I_inc(I_inc; tp = "Taiwan")
		append!(Taiwan_ind, flag_Taiwan)

		flag_China = check_cond_single_I_inc(I_inc; tp = "China")
		append!(China_ind, flag_China)
	end
	JLD2.jldsave(vars.path_rec,
		r_obj = r_obj, Korea_ind = Korea_ind, Taiwan_ind = Taiwan_ind, China_ind = China_ind)
end

function calc_imp_prob(r_obj::NamedTuple; ind_num = [])
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

function add_hatched_pattern!(pl::Plots.Plot, prob::Vector)
	ind = argmax(prob .== 1.0)
	x = [0, prob[ind] * 100]
	y = [ind, ind] .- 0.40
	plot!(pl, x, y, fillrange = y .+ 0.85,
		fillstyle = :/, label = false, color = :black,
	)
end

function plot_imp(df::DataFrame, col_imp, col_label;
	title = "", yticks = true, kwds...)::Plots.Plot
	n_cnt = nrow(df)
	yticks_label = yticks == true ?
				   (1:nrow(df), df[:, col_label]) :
				   (1:N_CNT_ECL_JPN, ["" for i in 1:N_CNT_ECL_JPN])

	pl = bar(df[:, col_imp] * 100,
		yticks = yticks_label,
		title = title;
		kwds...,
	)
	add_hatched_pattern!(pl, df[:, col_imp])
	return pl
end

function insert_imp_region_legend!(pl::Plots.Plot; legend=(0.3, 0.2))
	cmap = palette(:default)
	blue, green, red = cmap[1], cmap[3], cmap[7]
	bar!(pl, [0 0 0],
		legend = legend,
		labels = [" Eastern Asia" " South-eastern Asia" " Other Asia"],
		color = [blue green red],
		lw = 1,
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		orientation = :horizontal,
	)
end

function add_info_to_df_for_barplot(df_mer::DataFrame; sort_col = :prob1)
	df_mer = insert_obs_date(df_mer)
	insert_colormap_info!(df_mer)
	insert_tick_labels!(df_mer)
	df_mer = @orderby(df_mer, $sort_col)
	return df_mer
end

function plot_imp_prob_global_fs_num_cnts(
	vars::Vector{VisVars}, titles::Vector{String};
	sort_col = :prob1, kwds = nothing,
)
	kwds = deepcopy(kwds)
	col_names = [:prob1, :prob2, :prob3, :prob4]
	# Prepare merged dataframes
	I_incs = [var.I_inc for var in vars]
	println(size.(I_incs))
	df_imps = filtered_imp_prob.(I_incs, col_names)
	df_mer = create_imp_prob_dataframe(df_imps...; sort_col = sort_col)
	df_mer = add_info_to_df_for_barplot(df_mer; sort_col = sort_col)

	kwds[:color] = colors = df_mer[:, :color]

	yticks = [true, false, false, false]
	pls = broadcast((col, title, yticks) ->
			plot_imp(df_mer, col, :tick_label; title = title, yticks = yticks, kwds...),
		col_names, titles, yticks,
	)
	insert_imp_region_legend!(pls[1])

	# global final size epidemics,
	ylabels = ["Proportion (%)", "", "", ""]
	pls_fs = broadcast((I_inc, ylabel) -> global_final_size(I_inc; ylabel = ylabel),
		I_incs, ylabels)
	pls_imp = broadcast((I_inc, ylabel) ->
			barplot_I_cum_countries(I_inc, cutoff = 0, ylabel = ylabel),
		I_incs, ylabels)
	# Imp countries
	annotate!(pls[1], (-0.7, 1.03), text("A", :black, :left, 20))
	annotate!(pls_fs[1], (-0.7, 0.9), text("B", :black, :left, 20))
	annotate!(pls_imp[1], (-0.7, 0.9), text("C", :black, :left, 20))

	# adjust figures
	pl = [pls..., pls_fs..., pls_imp...]
	l = @layout [
		a{0.8h} b c d;
		e f g h;
		i j k l
	]
	plot(pl..., layout = l, size = (900, 1000), fmt = :png, dpi = 300)
end

"""
...
Args
- `path_rec`: Recording path for summarised statistics.
...
"""
function plot_imp_prob_global_fs_num_cnts_large_data(
	vars::VisVars, titles;
	sort_col = :prob1, kwds = nothing,
)
	sum_obj = load(vars.path_rec)
	r_obj = sum_obj["r_obj"]

	Korea_ind = [i for (i, x) in enumerate(sum_obj["Korea_ind"])
				 if x == true]
	Taiwan_ind = [i for (i, x) in enumerate(sum_obj["Taiwan_ind"])
				  if x == true]
	China_ind = [i for (i, x) in enumerate(sum_obj["China_ind"])
				 if x == true]

	# Importation prob.
	df_mer = @chain CSV.read(PATH_POP, DataFrame) begin
		@transform begin
			:prob1 = calc_imp_prob(r_obj)
			:prob2 = calc_imp_prob(r_obj; ind_num = Korea_ind)
			:prob3 = calc_imp_prob(r_obj; ind_num = Taiwan_ind)
			:prob4 = calc_imp_prob(r_obj; ind_num = China_ind)
		end
		@subset :iso_code .!= "JPN"
		create_imp_prob_dataframe(_, sort_col = sort_col)
	end

	# Edit dataframe for visualisation.
	df_mer = add_info_to_df_for_barplot(df_mer; sort_col = sort_col)
	kwds[:color] = colors = df_mer[:, :color]

	col_names = [:prob1, :prob2, :prob3, :prob4]
	yticks = [true, false, false, false]
	pls = broadcast((col, title, yticks) ->
			plot_imp(df_mer, col, :tick_label; title = title, yticks = yticks, kwds...),
		col_names, titles, yticks,
	)
	insert_imp_region_legend!(pls[1])

	# global final size epidemics,
	fs = r_obj.fs .|> Float64

	ylabels = ["Proportion (%)", "", "", ""]
	conds = [isreal.(fs), Korea_ind, Taiwan_ind, China_ind]
	pls_fs = broadcast((ind, ylabel) -> global_final_size(fs[ind]; ylabel = ylabel),
		conds, ylabels)
	pls_imp = broadcast((ind, ylabel) ->
			barplot_I_cum_countries(r_obj.num_cnt[ind], ylabel = ylabel),
		conds, ylabels)

	annotate!(pls[1], (-0.7, 1.03), text("A", :black, :left, 20))
	annotate!(pls_fs[1], (-0.7, 0.9), text("B", :black, :left, 20))
	annotate!(pls_imp[1], (-0.7, 0.9), text("C", :black, :left, 20))

	# adjust figures
	pl = [pls..., pls_fs..., pls_imp...]
	l = @layout [
		a{0.8h} b c d;
		e f g h;
		i j k l
	]
	plot(pl..., layout = l, size = (900, 1000), fmt = :png, dpi = 300)
end

function print_number_of_valid_sims(vars::VisVars)
	sum_obj = load(vars.path_rec)
	r_obj = sum_obj["r_obj"]
	println("# of all sims: ", length(r_obj.fs))
	println("# of Korea: ", sum(sum_obj["Korea_ind"]))
	println("# of Taiwan: ", sum(sum_obj["Taiwan_ind"]))
	println("# of China: ", sum(sum_obj["China_ind"]))
end

function quantiles_over_week(mat::Matrix, p::Vector{Float64})::Dict
	q_dict = Dict()
	for p_ in p
		qs = []
		for i in 1:size(mat)[2]
			q = quantile(mat[:, i], p_)
			push!(qs, q)
		end
		q_dict[p_] = qs
	end
	return q_dict
end

function quantile_importation_dates(I_inc)
	# TODO: Fetch only country dict
	INTER_SIM_BASE = return_inter_sim_base()
	country_dict = INTER_SIM_BASE.country_dict
	df_quantile = DataFrame(Dict(:country => [], :q05 => [], :q50 => [], :q95 => []))
	n_sim, n_days, n_cnt = size(I_inc)

	for ind_cnt in 1:n_cnt
		I_inc_cnt = I_inc[:, :, ind_cnt]
		flag_mat = I_inc_cnt[:, :] .> 0
		initial_dates = [findfirst(flag_mat[i, :]) for i in 1:n_sim]
		initial_dates = filter(x -> x != nothing, initial_dates)
		if length(initial_dates) == 0
			push!(df_quantile, (country_dict[ind_cnt], NaN, NaN, NaN))
			continue
		end
		qs = @pipe quantile.(Ref(initial_dates), [0.05, 0.50, 0.95]) .|> trunc(Int, _)
		push!(df_quantile, (country_dict[ind_cnt], qs...))
	end
	return df_quantile
end

"""This function uses I_inc in each VisVars.
"""
function plot_1st_importation_scenarios(var_lis::Vector{VisVars}, titles)
	I_incs = [var.I_inc for var in var_lis]
	df_imps = filtered_imp_prob.(I_incs, [:prob1, :prob2, :prob3, :prob4])
	df_mer = create_imp_prob_dataframe(df_imps...; sort_col = :prob1)
	df_mer = insert_obs_date(df_mer)
	sort!(df_mer, :prob1)

	imp_qs = quantile_importation_dates.(I_incs)
	yticks = [true, false, true, false]
	label = [true, false, false, false]
	pls = broadcast(
		(imp_qs, title, ytick, label) ->
			imp_date_vis(imp_qs, df_mer; title = title, yticks = ytick, label = label),
		imp_qs, titles, yticks, label,
	)

	plot!(pls[1], left_margin = 15Plots.mm, right_margin = 2Plots.mm)
	plot!(pls[2], left_margin = 0Plots.mm, right_margin = 5Plots.mm)
	plot(pls...,
		layout = @layout[a b; c d], size = (1000, 1500),
		fmt = :png,
	)
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

function imp_date_vis(imp_dates::DataFrame, df_mer::DataFrame; yticks = true, title = "",
	label = true, sort_col = :prob1,
)
	REPORT_DELAY_DAYS = 12
	DATE_NULL = Date(2025, 12, 31)

	cols_qs = [:q05, :q50, :q95]
	df_imp = copy(imp_dates)
	for col in cols_qs
		df_imp[!, col] = coalesce.(df_imp[!, col], 1000)
		# For reporting delay including 3 days incubation period and
		# 9 days onset to medical attendance.
		df_imp[!, col] .+= REPORT_DELAY_DAYS
	end
	# NOTE: JPN has dropped here.
	df_imp = leftjoin(df_mer, df_imp, on = :location_old => :country)
	insert_tick_labels!(df_imp)
	sort!(df_imp, [sort_col])

	x_min, st_date, x_max = Date(2022, 12, 25), Date(2023, 1, 16), Date(2024, 12, 31)
	xticks = [st_date + Month(i - 1) for i in 1:24]
	yticks = yticks == true ? (1:N_CNT_ECL_JPN, df_imp[:, :tick_label]) : (1:N_CNT_ECL_JPN, "")
	pl = plot(
		xlim = [x_min.instant.periods.value,
			x_max.instant.periods.value],
		ylim = [0, N_CNT_ECL_JPN + 1],
		xrotation = 90,
		xticks = (xticks, xticks),
		yticks = yticks,
		xlabel = "First importation date since 2023",
		title = title,
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		size = (700, 900),
	)

	label1 = "With ≥1 importations"
	label1_non = "Without importations"
	label2 = "Observed"
	if label == false
		label1 = label1_non = label2 = :none
	end
	blue, orange, green = palette(:default)[[1, 2, 3]]

	# Draw credible intervals.
	ys = 1:N_CNT_ECL_JPN
	dy = 0.25
	qs_lmh_date = [st_date .+ Day.(qs)
				   for qs in eachcol(df_imp[:, cols_qs])]
	x_l, x_med, x_h = reshape.(qs_lmh_date, 1, :)

	label1_lis = vcat(["", "", "", "", ""], label1, label1_non, ["" for i in 8:N_CNT_ECL_JPN])
	cond = df_imp[:, :first_repo] .== DATE_NULL
	color = @pipe [f == true ? blue : :blue for f in cond] |> reshape(_, 1, :)
	plot!(pl, [x_l; x_h], [ys'; ys'],
		color = color, linealpha = 0.3, lw = 8,
		label = reshape(label1_lis, 1, :),
	)

	# Draw median marker
	dy = 0.1
	plot!(pl, [x_med; x_med], [(ys .- dy)'; (ys .+ dy)'],
		color = color, linealpha = 0.4, lw = 3.5, label = "")
	# Draw observed dates.
	scatter!(pl, df_imp[:, :first_repo], 1:N_CNT_ECL_JPN,
		label = label2, color = orange, marker = :circle,
		markerstrokewidth = 0.6, markersize = 4,
	)
	return pl
end

function summarise_imp_prob(I_inc)
	df_imp_prob = CSV.read(PATH_POP, DataFrame)
	m_return, N0_MSM, country_dict = read_inter_country_data(PATH_FLIGHT, PATH_POP)
	# get imp prob
	n_sim, _, n_country = size(I_inc)
	cnt = @pipe sum(I_inc, dims = 2)[:, 1, :] .|>
				(x -> x > 0) |>
				sum(_, dims = 1) |>
				(x -> x[1, :])
	prob = cnt / n_sim

	df_imp_prob[!, :MSM_pop] .= N0_MSM
	df_imp_prob[!, :imp_prob] .= prob
	return df_imp_prob
end

function save_inc_imp(vars::VisVars; nmax = 50_000)
	I_inc = get_metrix_from_each_sim(vars.path_sim_res, get_incidence_I; nmax = nmax)
	JLD2.save_object(vars.path_inc, I_inc)
	df_imp_prob = summarise_imp_prob(I_inc)
	CSV.write(vars.path_imp, df_imp_prob)
end

function save_filtered_imp(vars::VisVars; nmax = 50_000)
	load_I_inc!(vars; check = true)
	I_inc = filter_I_inc(vars.I_inc; tp = "Korea", cut_off = 10)
	df_imp_prob = summarise_imp_prob(I_inc)
	CSV.write(vars.path_imp_fil, df_imp_prob)
end

function gen_based_Sankey_diagram(vars::VisVars, country_dict; nmax = 10000)
	path_objs = fetch_sim_paths(vars.path_sim_res)
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
	CSV.write(vars.path_exp_imp_gen, df_sum)
	return df_sum
end

"""
...
Args
- `n_sim`: Actual number of simulations.
	After conditioned, the file number is not matched with simulation number.
...
"""
function conditional_beta_SAR(vars::VisVars, cut_off, tp; n_sim = 50_000)
	load_I_inc!(vars; check = true)
	I_inc = vars.I_inc
	cond = filter_I_inc_cond(I_inc; cut_off = cut_off, tp = tp)

	# Prepare the simulated beta (whether or not values were saved or not.)
	path_objs = fetch_sim_paths(vars.path_sim_res)
	βs = ThreadsX.collect(JLD2.load(path, "β") for path in path_objs)
	println("Number of valid samples: ", length(βs))

	βs_qs = quantile(βs, [0.025, 0.5, 0.975])
	βs_r = round.(βs_qs, digits = 2)
	println("β : $(βs_r[2]) ($(βs_r[1]), $(βs_r[3]))")

	SAR_qs = @pipe (1 .- exp.(-βs_qs)) .* 100 .|> round(_, digits = 3)
	println("SAR : $(SAR_qs[2]) ($(SAR_qs[1]), $(SAR_qs[3]))")
end

function weekly_incidence_figure(jpn_weekly::Matrix, title;
	xlabel = "", ylabel = "", vspan=[-10, 23],
	ylim = [0.9, 30], yticks, xlim,
)
	ps = [0.25, 0.50, 0.75, 0.95]
	#jpn_weekly = one_country_I_inc_to_weekly(I_inc[:, :, ind0_cnt])
	q_dict = quantiles_over_week(jpn_weekly, ps)
	n_week = size(jpn_weekly)[2]
	pl = plot(
		yaxis = :log10,
		lw = 4,
		xlim = xlim,
		ylim = ylim,
		yticks = yticks,
		titlefontsize = 14,
		legendfontsize = 10,
		legendtitlefontsize = 10,
		xtickfontsize = 10,
		ytickfontsize = 10,
		title = title,
		xlabel = xlabel, ylabel = ylabel,
		fmt = :png, dpi = 200,
		foreground_color_legend = nothing,
		background_color_legend = nothing,
	)
	vspan!(pl, vspan, color = :black, alpha = 0.2, label = :none)
	for p in ps
		lw = p == 0.5 ? 2 : 1
		ls = p == 0.95 ? :dot : :solid

		q = q_dict[p]
		p_lab = Int64(p * 100)
		plot!(pl, 1:n_week, q .+ 1,
			lw = lw, ls = ls, color = "blue",
			label = "$(p_lab)th",
			legend_title = "Quantile",
		)
	end
	pl
end

function plot_multiple_incidence_curve(vars::VisVars;
	ylim=[0.9, 2000], anno_pos_x=10, anno_pos_y = 500, exclude_obs_panel=false,
	vspan=[-10, 23]
)
	# Load data
	sum_obj = JLD2.load(vars.path_rec)
	r_obj = sum_obj["r_obj"]
	# Convert data to appropriate data type.
	jpn_weekly_all = mapreduce(permutedims, vcat, r_obj.jpn_weekly)
	Korea_ind = [i for (i, x) in enumerate(sum_obj["Korea_ind"]) if x == true]
	jpn_weekly_Korea = jpn_weekly_all[Korea_ind, :]


	kwds = (xlim = [0, 159], ylim = ylim,
		ylabel = "Weekly cases + 1",
		yticks = ([1, 10, 100, 1000], ["1", "10", "100", "1000"]))
	pl1 = weekly_incidence_figure(jpn_weekly_Korea, ""; vspan=vspan, kwds...)
	pl2 = weekly_incidence_figure(jpn_weekly_all, "";
		xlabel = "Epidemiological week from 16 January 2023\n(medical attendance date)",
		vspan=vspan,
		kwds...,
	)
	# "-2" is to be matched with the fitting period.
	annotate!(pl1, anno_pos_x, anno_pos_y, text("A", :black, :left, 24))
	annotate!(pl2, anno_pos_x, anno_pos_y, text("B", :black, :left, 24))

	# Panel C info
	if exclude_obs_panel == false
		# Prepare the mpox incidence in Japan.
		df_jpn = @subset CSV.read(PATH_OWID, DataFrame) :location .== "Japan"
		df_jpn_week = daily_cases_to_weekly(df_jpn, :date, :new_cases)
		@subset!(df_jpn_week, :year .>= 2023)
		@transform!(df_jpn_week, :epi_week = 1:nrow(df_jpn_week))
		println("Epidemiological week is up to $(nrow(df_jpn_week))")

		pl3 = bar(df_jpn_week[:, :epi_week], df_jpn_week[:, :weekly_cases] .+ 1;
			bar_width = 1.0, label = "",
			xlabel = "Epidemiological week from 1 January 2023",
			yaxis = :log10,
			xtickfontsize = 10, ytickfontsize = 10,
			kwds...,
		)
		vspan!(pl3, [-10, 23 - 2], color = :black, alpha = 0.2, label = :none)
		annotate!(pl3, 10, anno_pos_y, text("C", :black, :left, 24))
		layout = @layout [a; b; c]
		plot(pl1, pl2, pl3, layout = layout, size = (800, 800), dpi = 300)

	else
		plot!(pl2, left_margin=8Plots.mm, bottom_margin=3Plots.mm)
		layout = @layout [a; b]
		plot(pl1, pl2, layout=layout, size=(800, 500), dpi=300)
	end
end

function Plots.plot(var::VisVars; legend=(0.3, 0.2))
	load_I_inc!(var)
	df_mer = @pipe filtered_imp_prob(var.I_inc, :prob1) |>
				create_imp_prob_dataframe(_; sort_col = :prob1) |>
				add_info_to_df_for_barplot(_)
	imp_qs = quantile_importation_dates(var.I_inc)
	kwds = deepcopy(BAR_KWDS)
	kwds[:color] = df_mer[:, :color]
	kwds[:xlabel] = "Simulated importation probability (%)"

	pl_bar = plot_imp(df_mer, :prob1, :tick_label; kwds...)
	insert_imp_region_legend!(pl_bar; legend=legend)

	pl_space = plot(legend=false,grid=true, ticks=false, showaxis=true,
		foreground_color_subplot=:white)

	ylabel = "Proportion (%)"
	pl_fs = global_final_size(var.I_inc; ylabel=ylabel)
	pl_imp = barplot_I_cum_countries(var.I_inc, cutoff=0, ylabel=ylabel,
		xlabel="Number of countries where ≥1 importations", xrotation=0)

	plot!(pl_bar, topmargin=10Plots.mm)
	plot!(pl_imp, topmargin=10Plots.mm)
	annotate!(pl_bar, (-0.1, 1.03), text("A", :black, :left, 20))
	annotate!(pl_fs, (-0.1, 1.2), text("B", :black, :left, 20))
	annotate!(pl_imp, (-0.1, 1.2), text("C", :black, :left, 20))

	layout = @layout[a [b{0.1h}; c{0.4h};  d{0.4h}]]
	plot(pl_bar, pl_space, pl_fs, pl_imp, layout=layout, size=(900, 700))
end