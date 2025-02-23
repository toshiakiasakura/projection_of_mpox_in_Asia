# Model for a single country, containing the alive particle filter algorithm.

include("./distributions.jl")

struct NetworkParams
	γ1::Float64
	γ2::Float64
	days::Int64
	N0::Int64
	β::Float64
	α::Float64
	κ::Float64
	ks::Vector{Float64}
	Pk::Vector{Float64}
	n_comp::Int64
	I0_init::Int64
	ind0_k::Int64
	θ_deno::Float64
	function NetworkParams(;
		γ1 = 0.19, γ2 = 0.19, days = 365,
		N0 = 10_000, β = 0.6, α = 0.10, κ = 0.77,
		n_comp = 50, I0_init = 1, ind0_k = 25,
		θ_deno = 0.0, ks, Pk,
	)
		new(γ1, γ2, days, N0, β, α, κ, ks, Pk, n_comp, I0_init, ind0_k, θ_deno)
	end
end

Base.@kwdef mutable struct NetworkSIRModel
	days::Int64
	n_comp::Int64
	N::Vector{Int64} = fill(0, n_comp)
	S::Matrix{Int64} = fill(0, days, n_comp)
	E::Matrix{Int64} = fill(0, days, n_comp)
	I::Matrix{Int64} = fill(0, days, n_comp)
	R::Matrix{Int64} = fill(0, days, n_comp)
	I_new::Matrix{Int64} = fill(0, days, n_comp)
end

Base.@kwdef struct ScenarioParams
	α::Float64 # Parameters for sexual partner distribution.
	κ::Float64
	ks::Vector{Float64}
	Pk::Vector{Float64}
	max_ind_k::Int64 # Highest sexual activity level present in Japan.
	γ2::Float64 # 1 / Infectious period
	file_prefix::String
	σ_β::Float64 # AFP-related params.
	mixin::Dict{Any, Any}
end

# For copy purpose
function ScenarioParams(
	sp::ScenarioParams;
	α           = sp.α,
	κ           = sp.κ,
	ks          = sp.ks,
	Pk          = sp.Pk,
	max_ind_k   = sp.max_ind_k,
	γ2          = sp.γ2,
	file_prefix = sp.file_prefix,
	σ_β         = sp.σ_β,
	mixin       = sp.mixin,
)
	ScenarioParams(
		α           = α,
		κ           = κ,
		ks          = ks,
		Pk          = Pk,
		max_ind_k   = max_ind_k,
		γ2          = γ2,
		file_prefix = file_prefix,
		σ_β         = σ_β,
		mixin       = mixin,
	)
end

SIM_DAY = 8 # Fixed length of days for the AFP algorithm

SCENARIO1 = ScenarioParams(
	α = NATSAL4W_PARMS[1], κ = NATSAL4W_PARMS[2],
	ks = KS_4W, Pk = PK_4W, max_ind_k = 42,
	γ2 = 1 / 10, file_prefix = "natsal4w_posterior_inf10",
	σ_β = 0.2,
	mixin = Dict("β_init" => 1.0),
)

SCENARIO2 = ScenarioParams(
	α = NATSAL1Y_PARMS[1], κ = NATSAL1Y_PARMS[2],
	ks = KS_1Y, Pk = PK_1Y, max_ind_k = 50,
	γ2 = 1 / 10, file_prefix = "natsal1y_posterior_inf10",
	σ_β = 0.05,
	mixin = Dict("β_init" => 0.1),
)

SCENARIO3 = ScenarioParams(
	α = NATSAL4W_PARMS[1], κ = NATSAL4W_PARMS[2],
	ks = KS_4W, Pk = PK_4W, max_ind_k = 42,
	γ2 = 1 / 7, file_prefix = "natsal4w_posterior_inf7",
	σ_β = 0.3,
	mixin = Dict("β_init" => 1.2),
)

SCENARIO4 = ScenarioParams(
	α = NATSAL4W_PARMS[1], κ = NATSAL4W_PARMS[2],
	ks = KS_4W, Pk = PK_4W, max_ind_k = 42,
	γ2 = 1 / 14, file_prefix = "natsal4w_posterior_inf14",
	σ_β = 0.3,
	mixin = Dict("β_init" => 0.8),
)

SCENARIO1_CUT1000 = ScenarioParams(
	α = SCENARIO1.α, κ = SCENARIO1.κ, ks = SCENARIO1.ks,
	Pk = PK_4W_CUT1000, max_ind_k = 38,
	γ2 = SCENARIO1.γ2, file_prefix = "natsal4w_posterior_inf10_cut1000",
	σ_β = 0.2,
	mixin = Dict("β_init" => 1.0),
)

SCENARIO1_D0415 = let
	targetdata = CSV.read(PATH_JPN_0415, DataFrame)[:, :count_sum]
	ScenarioParams(SCENARIO1;
		file_prefix = "natsal4w_posterior_inf10_d0415",
		mixin = Dict("β_init" => 1.4, "targetdata" => targetdata),
	)
end

SCENARIO1_ASSORT = ScenarioParams(SCENARIO1;
	file_prefix = "natsal4w_posterior_inf10_assort",
	# proportion of MSM population is reduced to be 69.4%.
	mixin = Dict("β_init" => 1.0, "red_msm_prop" => 0.694),
)

function set_values(model::NetworkSIRModel,
	ind::Int64,
	S::Vector{Int64},
	E::Vector{Int64},
	I::Vector{Int64},
	R::Vector{Int64},
	I_new::Vector{Int64},
)
	model.S[ind, :] = S
	model.E[ind, :] = E
	model.I[ind, :] = I
	model.R[ind, :] = R
	model.I_new[ind, :] = I_new
	return nothing
end

function get_values(model::NetworkSIRModel, ind::Int64)
	return model.S[ind, :], model.E[ind, :], model.I[ind, :], model.R[ind, :]
end

function Base.copy(model::NetworkSIRModel)
	model_new = NetworkSIRModel(;
		days = model.days, n_comp = model.n_comp,
		N = model.N, S = copy(model.S), E = copy(model.E),
		I = copy(model.I), R = copy(model.R), I_new = copy(model.I_new),
	)
	return model_new
end

"""
	alive_particle_filter(
		p_dict::NamedTuple,
		targetdata::Vector{Int64},
		nparticles::Int64,
		pmf::Vector{Float64},
		β::Float64,
		p0::Float64;
		memory_efficient::Bool=false,
		Nmax::Float64=Inf,
		ll_pre::Float64=0.0,
		u::Float64=0.0,
		ll_prior::Float64=0.0,
	)

# Arguments
- `p::NetworkParams`: Parameters.
- `model::NetworkSIRModel`: Initialised model.
- `nparticles::Int64`: Particles for filter.

# Keywords
- `memory_efficient::Bool=false`:
	- if true, the past trajectory is not stored, and days of simulations should be 8.
	- if false, the past trajectory is saved.

# Returns
- `Dict`: Dictionary contains the each trajectory and used parameters.

# Note
- In p_dict, `θ_deno` and `ind0_k` is overwritten in the alive particle filter.
- β is prioritised over the value set in p_dict.

# References
- McKinley (2020), Bayesian Anal. 15(3): 839-870, DOI: 10.1214/19-BA1174
"""
function alive_particle_filter(
	p_dict::NamedTuple,
	targetdata::Vector{Int64},
	nparticles::Int64,
	pmf::Vector{Float64},
	β::Float64,
	p0::Float64;
	memory_efficient::Bool = false,
	Nmax::Float64 = Inf,
	ll_pre::Float64 = 0.0,
	u::Float64 = 0.0,
	ll_prior::Float64 = 0.0,
)
	T = length(targetdata)
	nt_total = 0
	status = "Success"

	b_mat = design_matrix(targetdata, p0)
	ll = -Inf
	ll_init = T * log(nparticles) + ll_prior
	# N+1 match is requred from the theory.
	nt_rec = [nparticles + 1 for _ in 1:T]

	# Avoid overwriting the old values.
	p_new = p_dict
	@reset p_new.β = β
	if (memory_efficient == true) & (p_new.days != 8)
		error("days should be 8 for memory efficient version.")
	end

	params, model = initialise_model(p_new, pmf)
	x::Vector{NetworkSIRModel} = [model for _ in 1:(nparticles+1)]
	ps::Vector{NetworkParams} = [params for _ in 1:(nparticles+1)]
	time = @elapsed begin
		for t in 1:T
			nt = 0
			ps_tmp = ps
			x_tmp = x
			t_week = memory_efficient == true ? 1 : t

			for j in 1:(nparticles+1)
				δ = 0
				for _ in 1:1.5e6
					if t == 1
						params, model = initialise_model(p_new, pmf)
						update_initial_value(params, model)
						for t_day in 2:7
							run_sim_one_step(params, model, t_day)
						end
						# since 8 is mapped to 1 when updating
						if memory_efficient == true
							set_values(model, 8,
								model.S[7, :], model.E[7, :], model.I[7, :],
								model.R[7, :], model.I_new[7, :],
							)
						end
					else
						ind = rand(1:nparticles)
						params::NetworkParams = ps[ind]
						model::NetworkSIRModel = copy(x[ind])
						update_model_1week(params, model, t_week;
							memory_efficient = memory_efficient)
					end
					nt += 1
					nt_total += 1

					# Stop criteria
					if nt_total > Nmax
						status = "Over Nmax"
						break
					end
					nt_rec[t] = nt
					ll = ll_init - sum(log.(nt_rec .- 1))
					if exp(ll - ll_pre) < u
						status = "Early stop"
						break
					end

					if memory_efficient == true
						I_obs7 = convert_Inew_degree_to_weekly_obs(model.I_new, 1; index = 2)
					else
						I_obs7 = convert_Inew_degree_to_weekly_obs(model.I_new, t_week)
					end
					# 50 represents unrealistic number, so does not account.
					if I_obs7[t_week] > 50
						continue
					end

					# b_mact index 1 represents 0 value, so +1 is required.
					if b_mat[t, I_obs7[t_week]+1] == 1
						#print(params.ind0_k, ", ")
						ps_tmp[j] = params
						x_tmp[j] = copy(model)
						δ = 1
						break
					end
				end
				if status == "Over Nmax"
					break
				end
				if status == "Early stop"
					break
				end
				if δ == 0
					error("Can not obtain valid particles.")
				end
			end
			ps = ps_tmp
			x = x_tmp
			nt_rec[t] = nt
			ll = ll_init - sum(log.(nt_rec .- 1))
		end
	end  # end for @elapsed
	return Dict(
		"loglikelihood" => ll, "models" => x,
		"params" => ps, "nt" => nt_rec,
		"b_mat" => b_mat, "p0" => p0,
		"elapsed_time" => time, "status" => status,
		"nt_total" => nt_total,
	)
end

"""Update model for 1 week.

# Keywords
- `memory_efficient`: If true, the day 8 is moved to the day 1, and
	updated simulation results are stored in day 2 to day 8.

# Returns
- Nothing: model object is internally edited.
"""
function update_model_1week(
	params::NetworkParams,
	model::NetworkSIRModel,
	t_week::Int64;
	memory_efficient::Bool = false,
)
	if memory_efficient == true
		if t_week != 1
			error("t_week is not correctly set")
		end

		model.S[1, :] = model.S[8, :]
		model.E[1, :] = model.E[8, :]
		model.I[1, :] = model.I[8, :]
		model.R[1, :] = model.R[8, :]
		model.I_new[1, :] = model.I_new[8, :]
		for t_day in 2:8
			run_sim_one_step(params, model, t_day)
		end
	else
		st = (t_week - 1) * 7 + 1
		fin = t_week * 7
		for t_day in st:fin
			run_sim_one_step(params, model, t_day)
		end
	end
end

"""
	initialise_model(args...) -> (NetworkParams, NetworkSIRModel)

Initialised model is returned with calculated θ_deno.
Use `update_initial_value` and `run_sim_one_step`.

# Arguments
- `p::NamedTuple`: Initialised parameters.
- `pmf::Vector{Float64}`: Probablity mass function for degrees
	given one infectious contact.
	`get_probs_given_a_contact(p.α, p.κ; max_ind_k=max_ind_k)`
"""
function initialise_model(p::NamedTuple, pmf::Vector{Float64})
	@unpack n_comp, N0, ks, Pk, γ2 = p
	# Allocate N based on the Pk

	@reset p.ind0_k = wsample(1:50, pmf, 1)[1]

	N = N0 .* Pk .|> round .|> Int64
	θ_deno = (
		ks[i] / 365 * N[i]
		for i in 1:n_comp
	) |> sum
	@reset p.θ_deno = θ_deno

	params = NetworkParams(; p...)
	model = NetworkSIRModel(n_comp = 50, N = N, days = p.days)
	return params, model
end

"""
	update_initial_value(args...) -> (NetworParams, NetworkSIRModel)

Initialise model using pre-existing model.
This function only update the value at time 1.

"""
function update_initial_value(
	p::NetworkParams,
	model::NetworkSIRModel,
)
	@unpack N = model
	S0 = copy(N)
	S0[p.ind0_k] = S0[p.ind0_k] - p.I0_init
	E0 = fill(0, p.n_comp)
	I0 = fill(0, p.n_comp)
	I0[p.ind0_k] = p.I0_init
	R0 = fill(0, p.n_comp)
	I_new0 = fill(0, p.n_comp)
	I_new0[p.ind0_k] = p.I0_init
	set_values(model, 1, S0, E0, I0, R0, I_new0)
end

"""
	run_sim_one_step(args...) -> NetworkSIRModel

Run simulations for one step. This function is for the alive particle filter.
"""
function run_sim_one_step(
	p::NetworkParams,
	model::NetworkSIRModel,
	t::Int64,
)
	@unpack γ1, γ2, days, N0, β, n_comp, θ_deno, ks = p

	S, E, I, R = get_values(model, t - 1)
	if sum(E) + sum(I) == 0
		I_new = fill(0, n_comp)
		set_values(model, t, S, E, I, R, I_new)
		return model
	end

	# Calculate θ
	θ_nume = (
		maximum([0.0, ks[i] / (365 * γ2) - 1]) * γ2 * I[i]
		for i in 1:n_comp
	) |> sum
	θ = θ_nume / θ_deno

	# SEIR model
	p_inf = β .* ks ./ 365 .* θ

	new_inf = 0
	rec_E = 0
	rec_I = 0
	try
		new_inf = rand_binom.(S, 1 .- exp.(-p_inf))
		rec_E = rand_binom.(E, 1 - exp(-γ1))
		rec_I = rand_binom.(I, 1 - exp(-γ2))
	catch e
		println(p)
		println(t)
		println(S)
		println(E)
		println(I)
		println(new_inf)
		println(rec_E)
		println(rec_I)
		error(e)
	end
	newS = S .- new_inf
	newE = E .+ new_inf .- rec_E
	newI = I .+ rec_E .- rec_I
	newR = R .+ rec_I
	set_values(model, t, newS, newE, newI, newR, new_inf)
end

function custom_AFP_MH(
	params::NamedTuple, targetdata::Vector{Int},
	nparticles::Int64, pmf::Vector{Float64},
	p0::Float64;
	iterations = 100, scenario_parms::ScenarioParams,
)::Chains

	function loglikelihood(β::Float64; ll_pre = -Inf, u = 0.0, ll_prior = 0.0)
		res_dic = alive_particle_filter(
			params, targetdata, nparticles, pmf, β, p0;
			memory_efficient = true, Nmax = 5000.0 * nparticles,
			ll_pre = ll_pre, u = u, ll_prior = ll_prior,
		)
		ll = res_dic["status"] == "Success" ? res_dic["loglikelihood"] : -Inf
		#if res_dic["status"] != "Success"
		#    println(res_dic["status"], res_dic["SAR"])
		#end
		return (ll, res_dic["status"])
	end
	@unpack σ_β = scenario_parms

	#prior_dist = Beta(1.5, 1.5)
	prior_dist = Truncated(Cauchy(0.0, 2.5), 0.0, Inf)
	status = fill("", iterations)
	accept = fill(0, iterations)
	p_β = fill(Inf, iterations)
	lls = fill(0.0, iterations)

	p_β[1] = scenario_parms.mixin["β_init"] # rand(prior_dist)
	ll_tmp, status[1] = loglikelihood(p_β[1])
	lls[1] = ll_tmp + Distributions.logpdf(prior_dist, p_β[1])

	@showprogress for i in 2:iterations
		u = rand()
		# Generate a candidate.
		p_β[i] = rand(Normal(p_β[i-1], σ_β))
		ll_prior = Distributions.logpdf(prior_dist, p_β[i])

		if ll_prior == -Inf # outside range is abruptly discarded.
			p_β[i] = p_β[i-1]
			lls[i] = lls[i-1]
			status[i] = "Proposal is outside range"
			continue
		end

		# Return value includes the prior probability.
		lls[i], status[i] = loglikelihood(p_β[i]; ll_pre = lls[i-1], u = u, ll_prior = ll_prior)

		# Then jump to the NEW POSITION
		if (status[i] == "Success") & (u < exp(lls[i] - lls[i-1]))
			accept[i] = 1
			continue
		else  # Stick with the CURRENT POSITION
			p_β[i] = p_β[i-1]
			lls[i] = lls[i-1]
		end
	end
	return Chains(p_β, ["β"],
		info = (loglikelihood = lls, status = status,
			accept = accept, params = params, p0 = p0, pmf = pmf,
			targetdata = targetdata,
		),
	)
end

function fit_β(scenario_parms::ScenarioParams;
	days = :none, iterations = 20000, nparticles = 100,
)
	@unpack γ2, α, κ, ks, Pk, max_ind_k, file_prefix, σ_β = scenario_parms

	targetdata, N0_MSM = read_constant_values()
	targetdata = haskey(scenario_parms.mixin, "targetdata") ?
				 scenario_parms.mixin["targetdata"] : targetdata
	N0_MSM = haskey(scenario_parms.mixin, "red_msm_prop") ?
			 trunc(Int64, scenario_parms.mixin["red_msm_prop"] * N0_MSM) : N0_MSM
	n_point = size(targetdata)[1]

	# base params for the simulation.
	p0 = 0.5
	days = days == :none ? n_ponit * 7 : days

	params = (γ1 = 1 / 3, γ2 = γ2, I0_init = 1, ind0_k = typemin(Int),
		α = α, κ = κ, β = NaN, n_comp = 50,
		days = days, N0 = N0_MSM,
		θ_deno = NaN, ks = ks, Pk = Pk,
	)
	pmf = get_probs_given_a_contact(params.α, params.κ; max_ind_k = max_ind_k)

	Random.seed!(5)
	@time chn = custom_AFP_MH(params, targetdata, nparticles, pmf, p0;
		iterations = iterations, scenario_parms = scenario_parms,
	)

	# Summarise sampling information
	println("!!!Note that these estimated values are not thinned.!!!")
	display(chn)
	plot(chn, fmt = :png) |> display
	chn.info[:status] |> countmap |> println
	chn.info[:accept] |> mean |> println

	# Save result
	dict = Dict("chn" => chn, "params" => params, "pmf" => pmf, "p0" => p0)
	now_str = get_now_str()
	path = "../tmp_afp/$(file_prefix)_$(now_str).jls"
	println(path)
	serialize(path, dict)
	return path
end

function show_traceplot(path; n_burn = 2000, n_thin = 10)
	res_dic = deserialize(path)
	chn = res_dic["chn"]
	# n_burn = length(chn)*0.1 |> floor |> Int64
	chn = chn[(n_burn+1):n_thin:end]
	chn |> display
	qs = quantile(chn; q = [0.025, 0.5, 0.975]) |> DataFrame
	q_lmh = @pipe qs[1, 2:4] |> values
	q_lmh_r = round.(q_lmh, digits = 2)
	println("β (95% conf.):	 $(q_lmh_r[2]) ($(q_lmh_r[1]), $(q_lmh_r[3]))")
	SAR_lmh = @. round(1 - exp(-q_lmh), digits = 3) * 100
	println("SAR (95% conf.): $(SAR_lmh[2]) ($(SAR_lmh[1]), $(SAR_lmh[3]))")

	plot(chn, fmt = :png, left_margin = 15Plots.pt, figtype = :png, bottom_margin = 6Plots.mm) |> display
end
