
Base.@kwdef mutable struct InterNetworkParams
    γ1::Float64 = 1 / 3
    γ2::Float64
    days::Int64 = 365 * 3
    N0_cnt::Vector{Int64}
    β::Float64
    α::Float64
    κ::Float64
    n_comp::Int64 = 50
    I0_init::Int64 = 1
    ind0_cnt::Int64 = 15 # Index for Japan.
    ind0_k::Int64 = 25
    n_country::Int64
    m::Array{Float64}
    c_const::Float64 = 1.0
    C_gen::Array{Float64} = fill(Inf, n_country)
end

Base.@kwdef mutable struct InterSEIRModel
    days::Int64
    n_comp::Int64
    n_country::Int64
    ks::Vector{Float64}
    Pk::Vector{Float64}
    Nk_cnt::Array{Int64,2} = fill(-100, n_country, n_comp)
    S::Array{Int64,3} = fill(-100, days, n_country, n_comp)
    E::Array{Int64,3} = fill(-100, days, n_country, n_comp)
    I::Array{Int64,3} = fill(-100, days, n_country, n_comp)
    R::Array{Int64,3} = fill(-100, days, n_country, n_comp)
    I_new::Array{Int64,3} = fill(-100, days, n_country, n_comp)
    import_event::Dict = Dict(:time => [], :import_cntry => [], :export_cntry => [], :count => [])
end

Base.@kwdef struct ResultInterCountrySEIRModel
    days::Int64
    n_comp::Int64
    n_country::Int64
    ks::Vector{Float64}
    Pk::Vector{Float64}
    N_MSM::Matrix{Int64}
    #I::Array{Int64, 3}
    I_new::Array{Int64,3}
    import_event::Dict = Dict(:time => [], :import_cntry => [], :export_cntry => [], :count => [])
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
        model.R[ind, :, :]
    )
end

"""

This function updates `Nk_cnt`, `C_gen`.
"""
function initialise_model(
    p_inter::InterNetworkParams,
)::Tuple{InterNetworkParams,InterSEIRModel}
    @unpack N0_cnt, C_gen, α, κ, n_country, n_comp = p_inter

    ks, Pk = degree_and_probability_weibull(α, κ)
    # Set degree distributions.
    Nk_cnt = fill(-100, n_country, n_comp)
    for i in 1:n_country
        Nk_cnt[i, :] .= round.(Pk .* N0_cnt[i]; digits=0) .|> Int64
    end

    # Calculate a C_gen.
    for i in 1:n_country
        C_gen[i] = sum(ks ./ 365 .* Nk_cnt[i, :])
    end

    p_inter.C_gen = C_gen

    model = InterSEIRModel(
        days=p_inter.days,
        n_comp=n_comp,
        n_country=n_country,
        ks=ks, Pk=Pk, Nk_cnt=Nk_cnt
    )
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
    p::InterNetworkParams, new_inf, θ, λ_import, time
)
    @unpack n_comp, n_country, c_const = p

    new_inf = sum(new_inf, dims=2) # Sum over degrees, k.
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
    t::Int64
)
    @unpack γ1, γ2, β, n_country, C_gen, ind0_cnt, m, n_comp, c_const = p_inter
    @unpack ks = model

    S, E, I, R = get_values(model, t - 1)
    # If exposed and infectious individuals are not present
    # Pass the same value to the next step.
    if sum(E) + sum(I) == 0
        I_new = fill(0, n_country, n_comp)
        set_values!(model, t, S, E, I, R, I_new)
    end

    C_inf = fill(Inf, n_country)
    # Calculate a C_inf
    for i in 1:n_country
        C_inf[i] = sum(
            max.(0, ks ./ (365 * γ2) .- 1) .* γ2 .* I[i, :]
        )
    end

    θ = C_inf ./ C_gen
    η = fill(Inf, n_country, n_country)
    for i in 1:n_country
        for j in 1:n_country
            η[j, i] = C_inf[j] * m[j, i] / C_gen[i]
        end
    end

    # Hazard considering both.
    λ = fill(Inf, n_country, n_comp)
    # Hazard for importation.
    λ_import = fill(Inf, n_country, n_country)
    for i in 1:n_country
        for j in 1:n_country
            λ_import[i, j] = c_const * (m[i, j] * θ[j] + η[j, i])
        end
        λi = β .* ks ./ 365 .* (θ[i] + sum(λ_import[i, :]))
        λ[i, :] = λi
    end

    new_inf = rand_binom.(S, 1 .- exp.(-λ))
    rec_E = rand_binom.(E, 1 .- exp.(-γ1))
    rec_I = rand_binom.(I, 1 .- exp.(-γ2))

    newS = S .- new_inf
    newE = E .+ new_inf .- rec_E
    newI = I .+ rec_E .- rec_I
    newR = R .+ rec_I
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


"""
    read_inter_country_data(arg...) -> (Matrix, Vector, Dict)

Returns
- `Matrix`: The volume of returning travellers per day.
- `Vector`: Population size
- `Dict`: Corresponding index and country.
"""
function read_inter_country_data(path_flight, path_pop; D_travel=7)
    # rows are destination countries, columns are departure countries

    flight_matrix = CSV.read(path_flight, DataFrame)
    flight_matrix = flight_matrix[:, 3:end] |> Matrix
    flight_matrix = transpose(flight_matrix)
    n_country = size(flight_matrix)[1]

    # To set the diagnal elements to be 0
    for i in 1:n_country
        flight_matrix[i, i] = 0
    end
    N_size = CSV.read(path_pop, DataFrame)
    N_size = sort(N_size, :location)

    country_dict = Dict()
    for i in 1:n_country
        country_dict[i] = N_size[i, :location]
    end

    # Data is per 1,000.
    N0_pop = N_size[:, :pop2022] .* 1_000
    # To convert the population into MSM population.
    N0_MSM = N0_pop .* 0.01
    N0_MSM = N0_MSM .|> round .|> Int64
    m_return = copy(flight_matrix) .|> Float64

    for row in 1:size(m_return)[1]
        m_return[row, :] = m_return[row, :] .* D_travel ./ (N0_pop[row] .* 365)
    end
    return (m_return, N0_MSM, country_dict)
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
    n_sim_samples=10,
    nparticles=200,
    n_thin=10,
    n_burn=100,
    n_remain=50,
)::Nothing
    if isdir(path_dir) == true
        error("To continue, delete existing $(path_dir)")
    else
        mkdir(path_dir)
    end

    targetdata, N0_MSM = read_constant_values()
    # Extract required information from APF MH fitting results.
    APF_MH_res = deserialize(path_APF_MH_res)
    chn = APF_MH_res["chn"]
    #n_burn = length(chn) * 0.1 |> floor |> Int64
    chn = chn[(n_burn+1):n_thin:end]

    βs = chn["β"].data[:, 1]
    println("Length of samples: $(length(βs))")
    pmf = APF_MH_res["pmf"]
    @unpack γ2, α, κ, N0 = APF_MH_res["params"]
    p0 = APF_MH_res["p0"]
    # This line is no longer needed since all samples are used.
    # βs_sampled = sample(βs, n_sim_samples; replace=false)

    # Set required parameters
    n_point = length(targetdata)
    days = n_point * 7

    params = (γ1=1 / 3, γ2=γ2, I0_init=1,
        α=α, κ=κ, n_comp=50,
        days=days, N0=N0,
        # β is set in the next step.
        β=NaN,
        # Updated in `alive_partifle_filter`
        θ_deno=NaN, ind0_k=typemin(Int),
    )

    @showprogress for (i, β) in enumerate(βs)
        sim_jpn = Dict(
            "params" => [], "models" => [],
            "status" => [], "loglikelihood" => [],
        )
        # params, pmf, β should come from fitting results.
        APF_res = alive_particle_filter(
            params, targetdata, nparticles,
            pmf, β, p0, Nmax=nparticles * 5000.0
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
    N0_cnt::Vector{Int64},
    m_return::Matrix,
)
    n_country = length(N0_cnt)
    @unpack γ2, β, α, κ, ind0_k = p_jpn
    p_inter = InterNetworkParams(
        γ2=γ2, β=β, α=α, κ=κ,
        ind0_k=ind0_k,
        N0_cnt=N0_cnt,
        n_country=n_country, m=m_return,
    )
    return p_inter
end

function run_simulation(
    p_inter::InterNetworkParams,
    model_jpn::NetworkSIRModel;
    Korea_cond::Bool=false,
)
    # Run simulations
    p_inter, model = initialise_model(p_inter)
    set_initial_value!(p_inter, model, model_jpn)
    for t in 2:p_inter.days
        run_sim_one_step(p_inter, model, model_jpn, t)
        if (Korea_cond == true) & (t == 24 * 7) # 24 is the fitting period.
            # 24 is the fitting period, 18 denotes Korea,
            # 10 represents cutoff value.
            if sum(model.I_new[begin:24*7, 18, :]) < 10 # 18 is Korea
                break
            end
        end
    end
    res = ResultInterCountrySEIRModel(
        days=model.days, n_comp=model.n_comp,
        n_country=model.n_country,
        ks=model.ks, Pk=model.Pk,
        N_MSM=model.Nk_cnt, # I=model.I,
        I_new=model.I_new,
        import_event=model.import_event,
    )
    return res
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
    sim_jpn_dir::String; Korea_cond=false,
)
    # Read the flight matrix, MSM population and country information
    path_flight = "../data/flight/selected_flight_matrix.csv"
    path_pop = "../data/pop_size_edit.csv"
    m_return, N0_MSM, country_dict = read_inter_country_data(path_flight, path_pop)
    n_country = size(m_return)[1]

    now_str = get_now_str()
    path = "../tmp_results/$now_str"
    mkdir(path)
    println(path)

    sim_jpn_paths = fetch_sim_paths(sim_jpn_dir)
    sim_ind = 0
    @showprogress for sim_path in sim_jpn_paths
        sim_jpn = JLD2.load(sim_path)
        # Repeated access to the same trajectories.
        n_len = sim_jpn["params"] |> length
        for ind in 1:n_len
            sim_ind += 1
            p_jpn = sim_jpn["params"][ind]
            model_jpn = sim_jpn["models"][ind]

            p_inter = map_jpn_to_inter_params(p_jpn, N0_MSM, m_return)
            res = run_simulation(p_inter, model_jpn; Korea_cond=Korea_cond)
            if (Korea_cond == true)
                # 24 is the fitting period, 18 denotes Korea,
                # 10 represents cutoff value.
                if sum(res.I_new[begin:24*7, 18, :]) < 10
                    continue
                end
            end
            JLD2.jldsave("$(path)/$(sim_ind).jld2", res=res, β=p_jpn.β)
        end
    end
end

function international_spread_trajectories(
    path_APF_MH_res;
    n_sim_samples=1000, nparticles=200
)
    sim_jpn = obtain_jpn_trajectories_for_each_β(
        path_APF_MH_res;
        n_sim_samples=n_sim_samples, nparticles=nparticles
    )
    run_and_save_intercountry_model(sim_jpn, n_sim_samples)
end

function get_metrix_from_each_sim(path::String, f; nmax::Int64=100_000)
    path_objs = fetch_sim_paths(path)
    n = length(path_objs)
    n = minimum([n, nmax])
    sim_res::Array{Any,1} = fill(Inf, n)
    @showprogress for i in 1:n
        res = JLD2.load(path_objs[i])["res"]
        sim_res[i] = f(res)
    end
    sim_mat = vec_matrix_to_matrix(sim_res)
    return sim_mat
end

function get_incidence_I(model::ResultInterCountrySEIRModel)
    return @pipe model.I_new |> sum(_, dims=3)[:, :, 1]
end

function compare_filenames(a::AbstractString, b::AbstractString)
    # Extract numerical part of filenames
    num_a = parse(Int, match(r"(\d+)", basename(a)).match)
    num_b = parse(Int, match(r"(\d+)", basename(b)).match)
    # Compare numerical parts
    return num_a < num_b
end

function fetch_sim_paths(path::String)
    paths = glob("$(path)/*.jld2")
    filter!(x -> (x ≠ "$(path)/sim_jpn.jld2"), paths)
    sorted_filenames = sort(paths, lt=compare_filenames)
    return sorted_filenames
end

function summarise_imp_prob(I_inc)
    df_imp_prob = CSV.read(path_pop, DataFrame)
    # get imp prob
    n_sim, _, n_country = size(I_inc)
    cnt = @pipe sum(I_inc, dims=2)[:, 1, :] .|>
                (x -> x > 0) |>
                sum(_, dims=1) |>
                (x -> x[1, :])
    prob = cnt / n_sim

    df_imp_prob[!, :MSM_pop] .= N0_MSM
    df_imp_prob[!, :imp_prob] .= prob
    return df_imp_prob
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
        q05 = @pipe quantile(initial_dates, 0.05) |> trunc(Int, _)
        q50 = @pipe quantile(initial_dates, 0.50) |> trunc(Int, _)
        q95 = @pipe quantile(initial_dates, 0.95) |> trunc(Int, _)
        push!(df_quantile, (country_dict[ind_cnt], q05, q50, q95))
    end
    return df_quantile
end
