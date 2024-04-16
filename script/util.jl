using Accessors
using CategoricalArrays
using CSV
using DataFrames
using Dates
using Distributions
using Glob
using JLD2
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
using TimeSeries
using QuadGK
using WeibullParetoDist
using XLSX

get_now_str() = @pipe now() |> Dates.format(_, "yyyymmdd_HHMMSS")

function weekly_num_cases(tab::DataFrame, col::String)
    tab_time = TimeArray(tab, timestamp=Symbol(col))
    t = timestamp(tab_time)
    st = t[1]
    fin = t[end]
    tab_vis = DataFrame(day=[], count_sum=[])
    for i in 0:(fin-st).value
        d_next = st + Day(i)
        if d_next in t
            c = values(tab_time[d_next])[1]
        else
            c = 0
        end
        push!(tab_vis, [d_next, c])
    end
    tab_vis = TimeArray(tab_vis, timestamp=Symbol("day"))
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
function convert_Inew_degree_to_weekly_obs(I_new::Array{Int64,2}, n_point; index=1)
    I_obs = sum(I_new, dims=2)[:, 1]
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
    int = quadgk(target, wb.lower, 1_000_000, rtol=1e-8)[1]
    R = SAR * int / m
    return R
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
    r = innerjoin(r1, r2, on=:parameters)
    return r
end

function read_constant_values()
    path = "../tmp_fix_results/mpox_japan_weekly_data_0707.csv"
    targetdata = CSV.read(path, DataFrame)[:, :count_sum]
    N0_MSM = 1_239_517
    #N0_MSM_tokyo = 140_476
    #n_point = size(targetdata)[1]
    return targetdata, N0_MSM
end

function read_observed_imp_data()::DataFrame
    path = "../data/mpox_Asia_importation_date.xlsx"
    df_obs = @pipe XLSX.readtable(path, "Sheet2", "A:F"; first_row=1) |>
                   DataFrame .|>
                   coalesce(_, 0)
end