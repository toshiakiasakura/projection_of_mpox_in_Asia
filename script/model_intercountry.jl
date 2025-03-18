# Simulate international spread of mpox across 42 countires.

include("./data_path.jl")
include("./model.jl")

const IND_KOREA::Int64 = 18
const IND_JAPAN::Int64 = 15

Base.@kwdef mutable struct InterNetworkParams
	γ1::Float64 = 1 / 3
	γ2::Float64
	days::Int64 = 365 * 3
	N0_cnt::Vector{Int64}
	β::Float64
	α::Float64
	κ::Float64
	ks::Vector{Float64} # These two were modified part.
	Pk::Vector{Float64}
	n_country::Int64
	n_comp::Int64 = 50
	I0_init::Int64 = 1
	ind0_cnt::Int64 = IND_JAPAN # Index for Japan.
	ind0_k::Int64 = 25
	m::Matrix{Float64} # Flight matrix
	c_const::Float64 = 1.0
	C_gen::Array{Float64} = fill(Inf, n_country)
	fit_week_len::Int64 = 24
end

Base.@kwdef mutable struct InterSEIRModel
	days::Int64
	n_comp::Int64
	n_country::Int64
	Nk_cnt::Array{Int64, 2} = fill(-100, n_country, n_comp) # Degree specific population size.
	S::Array{Int64, 3} = fill(-100, days, n_country, n_comp)
	E::Array{Int64, 3} = fill(-100, days, n_country, n_comp)
	I::Array{Int64, 3} = fill(-100, days, n_country, n_comp)
	R::Array{Int64, 3} = fill(-100, days, n_country, n_comp)
	I_new::Array{Int64, 3} = fill(-100, days, n_country, n_comp)
	import_event::Dict = Dict(:time => [], :import_cntry => [], :export_cntry => [], :count => [])
end

function initialise!(model::InterSEIRModel;
	days::Int64, n_comp::Int64, n_country::Int64, Nk_cnt::Matrix)
	if model.days != days || model.n_comp != n_comp || model.n_country != n_country
		raise("Errors: Assumptions are wrong")
	end
	# Dimensions are the same; modify arrays in place
	model.Nk_cnt = Nk_cnt
	fill!(model.S, -100)
	fill!(model.E, -100)
	fill!(model.I, -100)
	fill!(model.R, -100)
	fill!(model.I_new, -100)

	# Reset the import_event dictionary
	for key in keys(model.import_event)
		model.import_event[key] = []
	end
end

Base.@kwdef mutable struct ResultInterCountrySEIRModel
	days::Int64
	n_comp::Int64
	n_country::Int64
	N_MSM::Matrix{Int64}
	#I::Array{Int64, 3}
	I_new::Matrix{Int64}
	import_event::Dict = Dict(:time => [], :import_cntry => [], :export_cntry => [], :count => [])
end

"""Data holder for running intercountry simulations
to avoid a repeated memory allocation.
"""
Base.@kwdef mutable struct SimBuffers
	n_country::Int64
	n_comp::Int64
	C_inf::Vector{Float64} = zeros(n_country)
	θ::Vector{Float64} = zeros(n_country)
	η::Matrix{Float64} = zeros(n_country, n_country)
	λ::Matrix{Float64} = zeros(n_country, n_comp) # Hazard considering both.
	λ_import::Matrix{Float64} = zeros(n_country, n_country) # Hazard for importation.
	new_inf::Matrix{Int64} = zeros(n_country, n_comp)
	rec_E::Matrix{Int64} = zeros(n_country, n_comp)
	rec_I::Matrix{Int64} = zeros(n_country, n_comp)
	newS::Matrix{Int64} = zeros(n_country, n_comp)
	newE::Matrix{Int64} = zeros(n_country, n_comp)
	newI::Matrix{Int64} = zeros(n_country, n_comp)
	newR::Matrix{Int64} = zeros(n_country, n_comp)
end

Base.@kwdef mutable struct InterSimParms
	m_return::Matrix{Float64}
	N0_MSM::Vector{Int64}
	n_country::Int64 = 42
	n_comp::Int64 = 50
	fit_week_len = 24 # Length of the Japanese incidence for fitting.
	country_dict::Dict = Dict()
	mixin::Dict = Dict()
	cutoff_inc::Int64 = 10
end

"""
	read_inter_country_data(arg...) -> (Matrix, Vector, Dict)

Returns
- `Matrix`: The volume of returning travellers per day.
- `Vector`: Population size
- `Dict`: Corresponding index and country.
"""
function read_inter_country_data(path_flight, path_pop; D_travel = 7)
	prop_MSM = 0.01

	# Rows are destination countries, columns are departure countries
	# in the original data and in the model,
	# rows are departure countires, columns are destination countries.
	flight_matrix = CSV.read(path_flight, DataFrame)
	flight_matrix = flight_matrix[:, 3:end] |> Matrix
	flight_matrix = transpose(flight_matrix)
	n_country = size(flight_matrix)[1]
	# To set the diagnal elements to be 0
	for i in 1:n_country
		flight_matrix[i, i] = 0
	end

	n_country = size(flight_matrix)[1]
	N_size = CSV.read(path_pop, DataFrame)
	N_size = sort(N_size, :location)

	country_dict = Dict()
	for i in 1:n_country
		country_dict[i] = N_size[i, :location]
	end

	# Data is per 1,000.
	N0_pop = N_size[:, :pop2022] .* 1_000
	# To convert the population into MSM population.
	N0_MSM = N0_pop .* prop_MSM
	N0_MSM = N0_MSM .|> round .|> Int64
	m_return = copy(flight_matrix) .|> Float64

	for row in 1:size(m_return)[1]
		@. m_return[row, :] = m_return[row, :] * D_travel / (N0_pop[row] * 365)
	end
	return (m_return, N0_MSM, country_dict)
end

function return_inter_sim_base()::InterSimParms
	m_return, N0_MSM, country_dict = read_inter_country_data(PATH_FLIGHT, PATH_POP)
	InterSimParms(
		m_return=m_return, N0_MSM=N0_MSM,
		country_dict=country_dict,
	)
end

function set_values!(model::InterSEIRModel,
	ind::Int64,
	S::Matrix{Int64},
	E::Matrix{Int64},
	I::Matrix{Int64},
	R::Matrix{Int64},
	I_new::Matrix{Int64},
)
	model.S[ind, :, :] = S
	model.E[ind, :, :] = E
	model.I[ind, :, :] = I
	model.R[ind, :, :] = R
	model.I_new[ind, :, :] = I_new
	return nothing
end

function get_values(model::InterSEIRModel, ind::Int64)
	return (
		model.S[ind, :, :],
		model.E[ind, :, :],
		model.I[ind, :, :],
		model.R[ind, :, :],
	)
end

"""

This function updates `Nk_cnt`, `C_gen`.
"""
function initialise_model(
	p_inter::InterNetworkParams;
	model::Union{Nothing, InterSEIRModel} = nothing,
)::Tuple{InterNetworkParams, InterSEIRModel}
	@unpack N0_cnt, C_gen, ks, Pk, n_country, n_comp = p_inter

	# Set degree distributions.
	Nk_cnt = fill(-100, n_country, n_comp)
	for i in 1:n_country
		Nk_cnt[i, :] .= round.(Pk .* N0_cnt[i]; digits = 0) .|> Int64
	end

	# Calculate a C_gen.
	for i in 1:n_country
		p_inter.C_gen[i] = sum(ks ./ 365 .* Nk_cnt[i, :])
	end

	if typeof(model) == InterSEIRModel
		initialise!(model; days = p_inter.days,
			n_comp = n_comp, n_country = n_country, Nk_cnt = Nk_cnt)
	else
		model = InterSEIRModel(
			days = p_inter.days,
			n_comp = n_comp,
			n_country = n_country,
			Nk_cnt = Nk_cnt,
		)
	end
	return p_inter, model
end

function set_initial_value!(
	p_inter::InterNetworkParams,
	model::InterSEIRModel,
	model_jpn::NetworkSIRModel,
)
	@unpack n_country, n_comp, ind0_cnt = p_inter
	S0 = copy(model.Nk_cnt)
	E0 = fill(0, n_country, n_comp)
	I0 = fill(0, n_country, n_comp)
	R0 = fill(0, n_country, n_comp)
	I_new0 = fill(0, n_country, n_comp)

	S0[ind0_cnt, :] = model_jpn.S[1, :]
	E0[ind0_cnt, :] = model_jpn.E[1, :]
	I0[ind0_cnt, :] = model_jpn.I[1, :]
	R0[ind0_cnt, :] = model_jpn.R[1, :]
	I_new0[ind0_cnt, :] = model_jpn.I_new[1, :]
	set_values!(model, 1, S0, E0, I0, R0, I_new0)
end

function set_import_events!(model::InterSEIRModel, time, impo, expo, count)
	append!(model.import_event[:time], time)
	append!(model.import_event[:import_cntry], impo)
	append!(model.import_event[:export_cntry], expo)
	append!(model.import_event[:count], count)
end

function record_import_events!(model::InterSEIRModel,
	p::InterNetworkParams, new_inf, θ, λ_import, time,
)
	@unpack n_comp, n_country, c_const = p

	new_inf = sum(new_inf, dims = 2) # Sum over degrees, k.
	for i in 1:n_country
		λi = sum(λ_import[i, :])
		p_import = λi / (θ[i] + λi)
		# If no probability of importation, skip.
		if (p_import == 0) | (θ[i] + λi == 0)
			continue
		end
		n_imp = rand_binom(new_inf[i], p_import)
		# If no importation is observed, skip.
		if n_imp == 0
			continue
		end

		λ_prob = λ_import[i, :] ./ λi
		n_imp_cnt = rand(Multinomial(n_imp, λ_prob), 1)[:, 1]
		for j in 1:n_country # export countries
			if n_imp_cnt[j] != 0
				set_import_events!(model, time, i, j, n_imp_cnt[j])
			end
		end
	end
end

"""
	run_sim_one_step(
		p_inter::InterNetworkParams,
		model::InterSEIRModel,
		model_jpn::NetworkSIRModel,
		t::Int64
	)
"""
function run_sim_one_step(
	p_inter::InterNetworkParams,
	model::InterSEIRModel,
	model_jpn::NetworkSIRModel,
	t::Int64,
	buffers::SimBuffers,
)
	@unpack γ1, γ2, β, n_country, C_gen, ind0_cnt, m, n_comp, c_const, ks = p_inter
	@unpack C_inf, θ, η, λ, λ_import, new_inf, rec_E, rec_I, newS, newE, newI, newR = buffers

	S, E, I, R = get_values(model, t - 1)
	# If exposed and infectious individuals are not present
	# Pass the same value to the next step.
	if sum(E) + sum(I) == 0
		fill!(new_inf, 0)
		set_values!(model, t, S, E, I, R, new_inf)
	end

	# Calculate a C_inf
	for i in 1:n_country
		C_inf[i] = sum(
            max.(0, ks ./ (365 * γ2) .- 1) .* γ2 .* I[i, :]
        )
	end

	@. θ = C_inf / C_gen
	for i in 1:n_country
		c_gen_i_rev = 1 / C_gen[i]
		for j in 1:n_country
			η[j, i] = C_inf[j] * m[j, i] * c_gen_i_rev
		end
	end

	# Hazard calculation
	for i in 1:n_country
		for j in 1:n_country
			λ_import[i, j] = c_const * (m[i, j] * θ[j] + η[j, i])
		end
		λ[i, :] = β .* ks ./ 365 .* (θ[i] + sum(λ_import[i, :]))
	end

	@. new_inf = rand_binom(S, 1 - exp(-λ))
	@. rec_E = rand_binom(E, 1 - exp(-γ1))
	@. rec_I = rand_binom(I, 1 - exp(-γ2))

	@. newS = S - new_inf
	@. newE = E + new_inf - rec_E
	@. newI = I + rec_E - rec_I
	@. newR = R + rec_I
	# Overwride the new value with the Japanese trajectories.
	if t <= model_jpn.days
		newS[ind0_cnt, :] = model_jpn.S[t, :]
		newE[ind0_cnt, :] = model_jpn.E[t, :]
		newI[ind0_cnt, :] = model_jpn.I[t, :]
		newR[ind0_cnt, :] = model_jpn.R[t, :]
		new_inf[ind0_cnt, :] = model_jpn.I_new[t, :]
	end
	record_import_events!(model, p_inter, new_inf, θ, λ_import, t)
	set_values!(model, t, newS, newE, newI, newR, new_inf)
end

function fetch_sim_paths(path::String)
	paths = glob("$(path)/*.jld2")
	filter!(x -> (x ≠ "$(path)/sim_jpn.jld2"), paths)
	sorted_filenames = sort(paths, lt = compare_filenames)
	return sorted_filenames
end

function compare_filenames(a::AbstractString, b::AbstractString)
	# Extract numerical part of filenames
	num_a = parse(Int, match(r"(\d+)", basename(a)).match)
	num_b = parse(Int, match(r"(\d+)", basename(b)).match)
	# Compare numerical parts
	return num_a < num_b
end

"""Obtain the Japanese transmission trajectories for each β value.
The whole sampled data is thinned in `n_thin`, and
from those thinned sampled, we randomly chose the next candidate values.

# Args:
- n_remain: Number of particles saved per one AFP simulation.

# Note:
From the current implementation, iterations are 12000,
so the number of thinned samples are 1200, and
among thme, 1000 samples are used for the subsequent simulations.
"""
function obtain_jpn_trajectories_for_each_β(
	path_APF_MH_res, path_dir;
	nparticles = 200,
	n_thin = 10,
	n_burn = 100,
	n_remain = 50,
	sc::Union{Nothing, ScenarioParams} = nothing,
)::Nothing
	if isdir(path_dir) == true
		error("To continue, delete existing $(path_dir)")
	else
		mkdir(path_dir)
	end
	targetdata, N0_MSM = read_constant_values()
	if sc != nothing
		targetdata = sc.mixin["targetdata"]
	end

	# Extract required information from APF MH fitting results.
	APF_MH_res = deserialize(path_APF_MH_res)
	chn = APF_MH_res["chn"]
	#n_burn = length(chn) * 0.1 |> floor |> Int64
	chn = chn[(n_burn+1):n_thin:end]

	βs = chn["β"].data[:, 1]
	println("Length of samples: $(length(βs))")
	pmf = APF_MH_res["pmf"]
	@unpack γ2, α, κ, N0, ks, Pk = APF_MH_res["params"]
	p0 = APF_MH_res["p0"]

	# Set required parameters
	n_point = length(targetdata)
	days = n_point * 7

	params = (γ1 = 1 / 3, γ2 = γ2, I0_init = 1,
		α = α, κ = κ, n_comp = 50,
		days = days, N0 = N0,
		β = NaN,
		# β is set in the next step.
		θ_deno = NaN, ind0_k = typemin(Int),
		ks = ks, Pk = Pk,
	)

	ThreadsX.foreach(enumerate(βs)) do (i, β)
		sim_jpn = Dict(
			"params" => [], "models" => [],
			"status" => [], "loglikelihood" => [],
		)
		# params, pmf, β should come from fitting results.
		APF_res = alive_particle_filter(
			params, targetdata, nparticles,
			pmf, β, p0, Nmax = nparticles * 5000.0,
		)
		if APF_res["status"] != "Success"
			error("Error! Increase nparticles.")
		end
		inds = sample(1:nparticles, n_remain)
		for ind in inds
			push!(sim_jpn["models"], APF_res["models"][ind])
			push!(sim_jpn["params"], APF_res["params"][ind])
			push!(sim_jpn["status"], APF_res["status"])
			push!(sim_jpn["loglikelihood"], APF_res["loglikelihood"])
		end
		JLD2.save("$(path_dir)/$(i).jld2", sim_jpn)
	end
end

function map_jpn_to_inter_params(
	p_jpn::NetworkParams,
	inter_parms::InterSimParms
)::InterNetworkParams
	@unpack γ2, β, α, κ, ks, Pk, ind0_k = p_jpn
	p_inter = InterNetworkParams(
		γ2 = γ2, β = β, α = α, κ = κ, ks = ks, Pk = Pk,
		ind0_k = ind0_k,
		N0_cnt = inter_parms.N0_MSM,
		n_country = inter_parms.n_country, m = inter_parms.m_return,
		fit_week_len = inter_parms.fit_week_len,
	)
	return p_inter
end

function run_simulation!(
	res::ResultInterCountrySEIRModel,
	p_inter::InterNetworkParams,
	model_jpn::NetworkSIRModel;
	Korea_cond::Bool = false,
	model_holder::Union{Nothing, InterSEIRModel} = nothing,
	buffers::Union{Nothing, SimBuffers} = nothing,
)::Nothing
	# Run simulations
	p_inter, model = initialise_model(p_inter; model = model_holder)
	set_initial_value!(p_inter, model, model_jpn)

	for t in 2:p_inter.days
		run_sim_one_step(p_inter, model, model_jpn, t, buffers)
		if (Korea_cond == true) & (t == p_inter.fit_week_len * 7) # 24 is the fitting period.
			# 24 is the fitting period, 18 denotes Korea,
			# 10 represents cutoff value.
			if sum(model.I_new[begin:(p_inter.fit_week_len*7), IND_KOREA, :]) < 10 # 18 is Korea
				break
			end
		end
	end

	# Record simulation results.
	if res.days != model.days || res.n_comp != model.n_comp || res.n_country != model.n_country
		raise("Errors: Assumptions are wrong.")
	end
	res.N_MSM = model.Nk_cnt
	res.I_new = sum(model.I_new, dims=3)[:, :, 1]
	for key in keys(model.import_event)
		res.import_event[key] = model.import_event[key]
	end

end

function prepare_initial_model_res_memory_storage(
	p_jpn::NetworkParams, inter_parms::InterSimParms)
	p_inter = map_jpn_to_inter_params(p_jpn, inter_parms)
	p_inter, model = initialise_model(p_inter)
	res = ResultInterCountrySEIRModel(
		days = model.days, n_comp = model.n_comp,
		n_country = model.n_country,
		N_MSM = model.Nk_cnt, # I=model.I,
		I_new = sum(model.I_new, dims=3)[:, :, 1],
		import_event = model.import_event,
	)
	return model, res
end

"""If the length of `sim_jpn` is smaller than `n_sim`,
repeatedly use the simulations.

...
Args:
- `Korea_cond::Bool`: If true,
	Early reject at 24 weeeks, and only save simulations
	when >= 10 incidence in Korea is observed.
...
"""
function run_and_save_intercountry_model(
	sim_jpn_dir::String, inter_parms::InterSimParms;
	Korea_cond = false,
)
	now_str = get_now_str()
	path = "../tmp_results/$now_str"
	mkdir(path); println(path)

	sim_jpn_paths = fetch_sim_paths(sim_jpn_dir)
	sim_ind = 0
	Threads.@threads for sim_path in sim_jpn_paths
		sim_jpn = JLD2.load(sim_path)

		# Prepare a struct instance to save memory.
		model_holder, res = prepare_initial_model_res_memory_storage(
			sim_jpn["params"][1], inter_parms)
		buffers = SimBuffers(
			n_country = inter_parms.n_country,
			n_comp = inter_parms.n_comp)

		# "params" contain 50 sets of trajectories.
		n_len = sim_jpn["params"] |> length
		# NOTE: This part is to increase the number of simulations
		# 	for forecasting 1st importation dates.
		#for _ in 1:10
			for ind in 1:n_len
				p_jpn = sim_jpn["params"][ind]
				model_jpn = sim_jpn["models"][ind]

				p_inter = map_jpn_to_inter_params(p_jpn, inter_parms)
				run_simulation!(res, p_inter, model_jpn;
					Korea_cond = Korea_cond,
					model_holder = model_holder,
					buffers = buffers)

				sim_ind += 1
				if (Korea_cond == true)
					if sum(res.I_new[begin:(inter_parms.fit_week_len*7), IND_KOREA, :]) < 10
						continue
					end
				end
				JLD2.jldsave("$(path)/$(sim_ind).jld2", res = res, β = p_jpn.β)
			end
		#end
	end
end
