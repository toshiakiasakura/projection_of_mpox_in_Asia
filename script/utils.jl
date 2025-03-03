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
