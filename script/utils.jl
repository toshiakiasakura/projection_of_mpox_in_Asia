using Accessors
using CategoricalArrays
using CSV
using DataFrames
using DataFramesMeta
using Dates
using Distributions
using FStrings
using Glob
using HypothesisTests
using JLD2
using LinearAlgebra
using MCMCChains
using Parameters
using Pipe
using Plots
using ProgressMeter
using Random
using Serialization
using SpecialFunctions
using StatsBase
using StatsPlots
using Test
using TimeSeries
using ThreadsX
using QuadGK
using WeibullParetoDist
using XLSX

include("./data_path.jl")
get_now_str() = @pipe now() |> Dates.format(_, "yyyymmdd_HHMMSS")

function weekly_num_cases(tab::DataFrame, col::String)
	tab_time = TimeArray(tab, timestamp = Symbol(col))
	t = timestamp(tab_time)
	st = t[1]
	fin = t[end]
	tab_vis = DataFrame(day = [], count_sum = [])
	for i in 0:(fin-st).value
		d_next = st + Day(i)
		if d_next in t
			c = values(tab_time[d_next])[1]
		else
			c = 0
		end
		push!(tab_vis, [d_next, c])
	end
	tab_vis = TimeArray(tab_vis, timestamp = Symbol("day"))
	tab_vis = collapse(tab_vis, week, first, sum)
	tab_vis = DataFrame(tab_vis)
	return tab_vis
end

"""Convert the matrix for I_new with degrees to weekly observed values.

# Arguments
- `I_new`: 2 dimentional matrix with degree.
- `n_point`: Number of weekly observational values. i.e. for 35 days, set 5.

# Keywords
- `index`: The starting date of summation.
"""
function convert_Inew_degree_to_weekly_obs(I_new::Array{Int64, 2}, n_point; index = 1)
	I_obs = sum(I_new, dims = 2)[:, 1]
	I_obs7 = []
	for i in 1:n_point
		st = index + (i - 1) * 7
		fin = index - 1 + i * 7
		I_obs7 = vcat(I_obs7, sum(I_obs[st:fin]))
	end
	return I_obs7
end

"""Convert the vector of matrix to 3 dimension matrix.
"""
function vec_matrix_to_matrix(vec)
	if isa(vec, Vector) == false
		return vec
	end
	n_vec = length(vec)
	n_r, n_c = size(vec[1])
	mat = fill(99999, length(vec), n_r, n_c)
	for i in 1:n_vec
		mat[i, :, :] = vec[i]
	end
	return mat
end

function R_eff_and_peak_relationship(I_cum::Float64,
	wb::LeftTruncatedWeibull,
	SAR::Float64)
	m = mean(wb)
	function target(x)
		S = pdf(wb, x) * exp(-I_cum * x / m)
		#nume = x*(x-1)*S
		nume = x * (x) * S
	end
	int = quadgk(target, wb.lower, 1_000_000, rtol = 1e-8)[1]
	R = SAR * int / m
	return R
end

"""Converts daily case data in a DataFrame to weekly case data.
"""
function daily_cases_to_weekly(df::DataFrame, date_col::Symbol, case_col::Symbol)::DataFrame
	if typeof(df[:, date_col]) != Vector{Date}
		df[!, date_col] = Date.(df[:, date_col], "dd/mm/yyyy")
	end
    df[!, :week] = Dates.week.(df[:, date_col])
    df[!,:year] = Dates.year.(df[:,date_col])

    # Group by year and week and sum the cases
    weekly_df = combine(groupby(df, [:year, :week]), case_col => sum => :weekly_cases)
    return weekly_df
end


"""
	ABC_threshold_for_each_value

# Returns
- `Int64`: Acceptable error given the probability, p0 based on Poisson dist.
"""
function ABC_threshold_for_each_value(xt::Int64, p0::Float64)
	vals = 0:50
	ϵ = 0
	ys = pdf.(Poisson.(vals), xt)
	for ϵ in 0:100
		x_range = (xt-ϵ):(xt+ϵ)
		x_range = [i for i in x_range if i >= 0]
		p_range = ys[x_range.+1] |> sum
		if p_range > p0
			return ϵ
		end
	end
	error("Not Found")
end

"""Return design matrix given accepted probabilities.

# Returns
- Matrix{Int64}:
	- First element represents time,
	- the second element's index - 1 is corresponding to the value.
	- Ex. b_mat[4,6] is 1, representing the time 4 and value 5 is accepted.
"""
function design_matrix(targetdata, p0)
	n_point = length(targetdata)
	b_mat = fill(0, n_point, 51)
	ϵs = ABC_threshold_for_each_value.(targetdata, p0)
	for i in 1:n_point
		x = targetdata[i]
		ϵ = ϵs[i]
		l = maximum([x - ϵ + 1, 1])
		inds = l:(x+ϵ+1)
		b_mat[i, inds] .= 1
	end
	return b_mat
end

function extract_chain_info_weibull_pareto(chn::Chains)::DataFrame
	r1 = summarize(chn, mean, median, var, std) |> DataFrame
	r2 = hpd(chn) |> DataFrame
	r = innerjoin(r1, r2, on = :parameters)
	return r
end

function read_constant_values()
	targetdata = CSV.read(PATH_JPN_0707, DataFrame)[:, :count_sum]
	N0_MSM = 1_239_517
	#N0_MSM_tokyo = 140_476
	#n_point = size(targetdata)[1]
	return targetdata, N0_MSM
end

function read_observed_imp_data()::DataFrame
	path = "../data/mpox_Asia_importation_date.xlsx"
	df_obs = @pipe XLSX.readtable(path, "Sheet2", "A:F"; first_row = 1) |>
				   DataFrame .|>
				   coalesce(_, 0)
end

"""Given dataframe, it tabulates to calculate the weekly numbers.
"""
function convert_to_weekly_case(df_::DataFrame, col::String)
	tab = @pipe groupby(df_, col) |> combine(_, :count => sum)
	sort!(tab, med_date)
	tab_7d = weekly_num_cases(tab, col)
	return tab_7d
end

Base.@kwdef mutable struct CategoricalDegree
	prop::Vector{Float64}
	bins::Vector{Real} # left-closed and right-open intervals.
	bins_label::Vector{String}
	zero_category::Bool
	period_y::Float64
	name::String
	ref::String
	sexuality::String
	participant::String
	partner_type::String
	mixin::Dict{Any, Any} = Dict()
end

KweeChoy2014 = let
	p = (28 + 29 + 10 + 50) / 273
	CategoricalDegree(
		prop = [1 - p, p],
		bins = [0, 11, Inf],
		bins_label = ["≤10", ">10"],
		zero_category = false,
		period_y = 0.5,
		name = "Kwee Choy 2014",
		ref = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4099220/pdf/IPID2014-236240.pdf",
		sexuality = "Homosexual and/or bisexual",
		participant = "Clients seeking VCT services",
		partner_type = "male or female",
		mixin = Dict(:location => "Vietnam", :sample_size => 273),
	)
end

# TODO: Remove this instant since this is no longer appropriate to be included..
JeffreyGrierson2013 = CategoricalDegree(
	prop = [3.2, 45.3, 14.9, 36.5] ./ 100, # Table 6
	bins = [0, 1, 6, 11, Inf],
	bins_label = ["0", "1-5", "6-10", "≥11"],
	zero_category = true,
	period_y = 1.0,
	name = "Jeffrey Grierson 2013",
	ref = "https://www.researchgate.net/publication/257749158_Networks_of_MSM_in_Indonesia_A_2-mode_study_of_MSM_and_sites_of_engagement",
	sexuality = "Homosexual and/or bisexual",
	participant = "From MSM network, convenient sampling",
	partner_type = "men only (Table 6)",
	mixin = Dict(:location => "Indonesia",
		:sample_size => 1329,
		:note => "Table 2 contains Bisexual proportion (40.3%),
			22.8% have a current relationship with a womane (Table 3)",
	),
)

# TODO: Remove this instant since this is no longer appropriate to be included..
AymanAssi2019 = CategoricalDegree(
	prop = [22.0, 58.0, 11.0, 9.0] ./ 100,
	bins = [0, 2, 6, 11, Inf],
	bins_label = ["0-1", "2-5", "6-10", "≥11"],
	zero_category = false,
	period_y = 0.25,
	name = "AmanAssi2019",
	ref = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6806001/",
	sexuality = "self-identifieid MSM",
	participant = "Attending sexual health centre",
	partner_type = "male or female",
	mixin = Dict(:location => "Lebanon", :sample_size => 2238),
)

# TODO: Delete this dataset since it is a convenient sampling with CBO.
HaochuLi2021 = CategoricalDegree(
	#prop = [30.3, 28.9, 34.8, 6.0] ./ 100,
	#bins = [0, 1, 2, 7, Inf],
	#bins_label = ["0", "1", "2-6", "≥7"],
	prop = [30.3 + 28.9, 34.8, 6.0] ./ 100,
	bins = [0.5, 2, 7, Inf],
	bins_label = ["≤1", "2-6", "≥7"],
	zero_category = false,
	period_y = 0.5,
	name = "Haochu Li 2021",
	ref = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7839183/",
	sexuality = "had sex with another male in the past 12 months",
	participant = "Convenient sampling with CBO",
	partner_type = "male_only",
	mixin = Dict(:location => "China",
		# Note: sample size is 419
		:sample_size => 419 - 127,  # excluding 0
		:sample_size0 => 121, # sample size for 1.
	),
)

# TODO: Remove this data
SarangJang2024 = CategoricalDegree(
	prop = [38.4, 49.6, 12.0] ./ 100,
	bins = [0.5, 2, 10, Inf],
	bins_label = ["≤1", "2-9", "≥10"],
	zero_category = false,
	period_y = 0.5,
	name = "Sarang Jang 2024",
	ref = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11359825/",
	sexuality = "self-reported engagement in oral or anal sex with at least one male partner in the previous year",
	participant = "online survey",
	partner_type = "male only",
	mixin = Dict(:location => "Republic of Korea",
		:sample_size => 1389, :sample_size0 => 491,
	),
)

# TODO: Remove this data
DapengZhang2007 = CategoricalDegree(
	prop = [17.5 + 24.7, 40.4, 17.4] ./ 100,
	bins = [0.5, 2, 6, Inf],
	bins_label = ["≤1", "2-5", "≥6"],
	zero_category = false,
	period_y = 0.5,
	name = "Dapeng Zhang 2007",
	ref = "https://pmc.ncbi.nlm.nih.gov/articles/PMC2598657/pdf/571.pdf",
	sexuality = "Men who were 18 years old and over and had oral or anal sex with males in the past year",
	participant = "online survey",
	partner_type = "male only",
	mixin = Dict(:location => "China (mainland)",
		:sample_size => 2011, :sample_size0 => 415 + 583,
	),
)

# TODO: Remove this data
SinHowLim2015 = CategoricalDegree(
	# TODO: After determining the visualisation, delete another.
	prop = normalize([29.3, 49.2, 21.5] ./ 100, 1),
	bins = [1, 2, 6, Inf],
	bins_label = ["1", "2-5", "≥6"],
	#prop = normalize([49.2, 21.5] ./ 100, 1),
	#bins = [2, 6, Inf],
	#bins_label = ["2-5", "≥6"],
	zero_category = false,
	period_y = 0.5,
	name = "ASIMM, Malaysia only", # "Sin How Lim 2015",
	ref = "https://pubmed.ncbi.nlm.nih.gov/25865907/",
	sexuality = "had sex with another male in the past 6 months",
	participant = "Web-based survey on gay community-related web tools",
	partner_type = "male only",
	mixin = Dict(:location => "Malaysia",
		:sample_size => 1235, :sample_size0 => 362,
		:study_name => "AIMSS study",
	),
)

# TODO: Remove this data
SinHowLim2012 = CategoricalDegree(
	# TODO: After determining the visualisation, delete another.
	prop = normalize([27.8, 47.0, 14.2, 9.4, 1.6] ./ 100, 1),
	bins = [1, 2, 6, 11, 51, Inf],
	bins_label = ["1", "2-5", "6-10", "11-50", "≥51"],
	#prop = normalize([47.0, 14.2, 9.4, 1.6] ./ 100, 1),
	#bins = [2, 6, 11, 51, Inf],
	#bins_label = ["2-5", "6-10", "11-50", "≥51"],
	zero_category = false,
	period_y = 0.5,
	name = "ASIMM", #"Sin How Lim 2012",
	ref = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4405155/",
	sexuality = "had sex with another male in the past 6 months",
	participant = "Web-based survey on gay community-related web tools",
	partner_type = "male only",
	mixin = Dict(:location => "Asian countries",
		:sample_size => 10_413, :sample_size0 => 2895,
		:study_name => "AIMSS study",
	),
)

function Plots.plot!(
	pl::Plots.Plot, cate_deg::CategoricalDegree;
	kwds...,
)
	@unpack prop, bins, bins_label, zero_category, name, period_y, mixin = cate_deg

	lower = bins[1]
	bins = deepcopy(bins)
	bins[1] = bins[1] == 0 ? 1 : bins[1]
	sample_size = mixin[:sample_size]
	if zero_category == true
		sample_size -= mixin[:sample_size0]
		prop = normalize(prop[2:end], 1)
		bins = bins[2:end]
		bins_label = bins_label[2:end]
	end

	α, κ = NATSAL4W_PARMS
	wb_4w = truncated(WeibullPareto(α, κ * (1 / period_y)^α); lower = lower)
	α, κ = NATSAL1Y_PARMS
	wb_1y = truncated(WeibullPareto(α, κ * (1 / period_y)^α); lower = lower)
	# TODO: Should be matched with half-year period in a certain way.
	# Less than one should be removed?
	p_pdf_4w = @pipe cdf.(wb_4w, bins) |> (x -> x[2:end] .- x[1:(end-1)])
	p_pdf_1y = @pipe cdf.(wb_1y, bins) |> (x -> x[2:end] .- x[1:(end-1)])
	normalize!(p_pdf_4w, 1)
	normalize!(p_pdf_1y, 1)
	@test (round(sum(prop), digits = 1) == 1.0)

	yerr_mat = create_yerr_from_prop_and_sample_size(prop, sample_size, 3)
	groupedbar!(pl, [prop p_pdf_4w p_pdf_1y] .* 100;
		xticks = (1:length(bins_label), bins_label),
		label = ["$(name), n=$(sample_size)" "Natsal 4-week" "Natsal 1-year"],
		yerr = yerr_mat, color = [1 2 3],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		kwds...,
	)
	annotate!(pl, (0.5, 1.1), text("$(mixin[:location])", :black, :center, 12))
	# TODO: Inset part, will be ignored.
	#plot!(pl, inset = bbox(0.5, 0.1, 0.35, 0.4), subplot = 2)
	#prop = normalize(prop[2:end], 1)
	#p_pdf_4w = normalize(p_pdf_4w[2:end], 1)
	#p_pdf_1y = normalize(p_pdf_1y[2:end], 1)
	#groupedbar!(pl, [prop p_pdf_4w p_pdf_1y] .* 100,
	#	inset = bbox(0.5, 0.1, 0.35, 0.4), subplot = 2,
	#)
	return pl
end

"""
Args:
- m: Number of groups for each x ticks.
"""
function create_yerr_from_prop_and_sample_size(
	prop::Vector, sample_size::Int64, m::Int64)::Matrix
	yerr_mat = fill((NaN, NaN), length(prop), m)
	value = [round(p * sample_size; digits = 0) for p in prop]
	for (i, p) in enumerate(prop)
		v = round(p * sample_size; digits = 0)
		yerr = confint(BinomialTest(v, sample_size))
		yerr = (p - yerr[1], yerr[2] - p) .* 100
		yerr_mat[i, 1] = yerr
	end

	return yerr_mat
end

function Plots.plot(cate_deg::CategoricalDegree; kwds...)
	pl = plot()
	plot!(pl, cate_deg; kwds...)
end

function Newman_assortativity_coef(mat::Matrix)::Float64
	e_sum = sum(mat * mat)
	r = (sum(diag(mat)) - e_sum) / (1 - e_sum)
	return r
end

function convert_50comp_to_EricChow2016_bins(mat_e50, ks_1y)
	println("----- Converting process info. -----")
	ks_3m = ks_1y ./ 4
	conds = [(ks_3m .< 2.0), (2.0 .<= ks_3m .< 4), (4 .<= ks_3m)]
	mat_e3 = zeros(3, 3)
	for i in 1:3, j in 1:3
		mat_e3[i, j] = sum(mat_e50[conds[i], conds[j]])
	end
	println("50*50 matrix : ", Newman_assortativity_coef(mat_e50))
	println("3*3 matrix : ", Newman_assortativity_coef(mat_e3))
	display(mat_e3 .* 100)
	println("Our Natsal 1Y model within 3 months, < 2 partner, 2-3 partners, >= 4 partners (%): ", round.(diag(mat_e3) .* 100, digits = 3))
	return mat_e3
end

function increase_low_assortativity(
	mat::Matrix, α1::Float64, α23::Float64 = 1.; norm=true)
	conds = [(1:11), (12:15), (16:33), (34:50)]
	mat_e50_α = copy(mat)
	mat_e50_α[conds[1], conds[1]] .*= α1
	mat_e50_α[conds[2], conds[2]] .*= α23
	if norm == true
		mat_e50_α = normalize(mat_e50_α, 1)
	end
	return mat_e50_α
end

function calculate_proportion_of_total_partnership(mat::Matrix)::Nothing
    C_norm = sum(mat)
    inds = [(1:11), (12:15), (16:50)]
    p_part1 = sum(mat[inds[1], inds[1]])/C_norm
    p_part2 = sum(mat[inds[2], inds[2]])/C_norm
    p_part3 = sum(mat[inds[3], inds[3]])/C_norm
    println(f"Percentage of total partnership for 1 to 1, 2-3 to 2-3, >4 to >4 pairs: {p_part1*100:.2f}%, {p_part2*100:.2f}%, {p_part3*100:.2f}%")
end

function calculate_proportion_of_individuals(mat::Matrix)::Nothing
    inds = [(1:11), (12:15), (16:50)]

    n_ind = sum(mat; dims=2) ./ KS_4W
    C_ind = sum(n_ind)
    p_inds = []
    for ind in inds
        mat_diag = mat[ind, ind]
        n_ind = (sum(mat_diag; dims=2) ./KS_4W[ind]) |> sumj
        push!(p_inds, n_ind/C_ind)
    end
    println(f"Percentage of inkakdividuals for 1 to 1 and 2-3 to 2-3, >4 to >4 pairs: {p_inds[1]*100:.2f}%, {p_inds[2]*100:.2f}%, {p_inds[3]*100:.2f}%")
end

function create_rewired_degree_distribution(Pk::Vector, α1::Float64, α2::Float64
	)::Tuple
	inds = [(1:11), (12:15), (16:50)]
	w1 = normalize(Pk[inds[1]], 1)
	w2 = normalize(Pk[inds[2]], 1)
	Pk_new = copy(PK_4W)
	Pk_new[inds[1]] .-= α1 .* w1
	Pk_new[inds[2]] .-= α2 .* w2
	@test isapprox(sum(Pk_new), 1 - α1 - α2; atol=0.01)
	normalize!(Pk_new, 1)
	return (Pk_new, w1, w2)
end

function create_total_and_prop_partnership_dist(Pk_new, w1, w2, ks, α1, α2)::Tuple
	Pk_new_r = normalize(Pk_new .* ks, 1)
	Pk_new_w1 = normalize(Pk_new_r[inds1], 1)
	Pk_new_w2 = normalize(Pk_new_r[inds2], 1)
	mat_e50_prop = (1 - α1 - α2) .* (Pk_new .* ks) .* reshape(Pk_new_r, 1, 50)
	mat_e50_homo1 = α1 .* (w1 .* ks[inds1]) .* reshape(Pk_new_w1, 1, length(w1))
	mat_e50_homo2 = α2 .* (w2 .* ks[inds2]) .* reshape(Pk_new_w2, 1, length(w2))
	mat_e50 = copy(mat_e50_prop)
	mat_e50[inds1, inds1] += mat_e50_homo1
	mat_e50[inds2, inds2] += mat_e50_homo2

	# Note: this normalisation to mat_e50_prop is only for visualisation purpose.
	normalize!(mat_e50_prop, 1)
	normalize!(mat_e50, 1)
	return (mat_e50_prop, mat_e50)
end
