# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Julia 1.8.3
#     language: julia
#     name: julia-1.8
# ---

include("./util.jl")
include("./vis_util.jl")
include("./model.jl")
include("./model_intercountry.jl")

# # Data preparation
# Note: Files in `../tmp_results/**` is created in 2j, and those file size is around 900GB for Natsal 4w IP10 (leaving all simulation resutls), and around 15GB for other scenarios (leaving only conditioned results).

β = 0.091
SAR = 1- exp(-β)

path_flight = "../data/flight/selected_flight_matrix.csv"
path_pop = "../data/pop_size_edit.csv"
m_return, N0_MSM, country_dict = read_inter_country_data(path_flight, path_pop)

# +
vars1_IP10 = VisVars(
    path_sim_res = "../tmp_results/20240411_003226", # 50_000, conditional only with β,
    suffix = "path1_IP10"
)

vars5_IP10 = VisVars(
    path_sim_res = "../tmp_results/20240411_005101_fil", # 50_000 with Korea condition.
    suffix = "path5_IP10"
)

vars5_IP7 = VisVars(
    path_sim_res = "../tmp_results/20240411_013445", # 50_000, only conditional one with β,
    suffix = "path5_IP7"
)

vars5_IP14 = VisVars(
    path_sim_res = "../tmp_results/20240411_005430", # 50_000, only conditional one with β,
    suffix = "path5_IP14"
)
# -

# ## Filter the file  

vars5_IP10_ori = VisVars(
    path_sim_res = "../tmp_results/20240411_005101", # 50_000 without any conditions.
    suffix = "path5_IP10"
)
path_rec = "../tmp_results/summary_quantity.jld2"

move_eligible_file_to_directory(vars5_IP10_ori, "../tmp_results/20240411_005101_fil")

record_summary_statistics(vars5_IP10_ori, path_rec)

record_summary_statistics(vars5_IP10_ori, "../tmp_results/summary_quantity_test.jld2"; nmax=50000)

# # Save Inc and imp. prob.

nmax=50_000

function save_inc_imp(vars::VisVars)
    I_inc = get_metrix_from_each_sim(vars.path_sim_res, get_incidence_I; nmax=nmax)
    JLD2.save_object(vars.path_inc, I_inc)
    df_imp_prob = summarise_imp_prob(I_inc)
    CSV.write(vars.path_imp, df_imp_prob)
end

save_inc_imp(vars5_IP10)

save_inc_imp(vars1_IP10)
save_inc_imp(vars5_IP7)
save_inc_imp(vars5_IP14)
save_inc_imp(vars5_IP10)

function save_filtered_imp(vars::VisVars)
    load_I_inc!(vars; check=true)
    I_inc = filter_I_inc(vars.I_inc; tp="Korea", cut_off=10)
    df_imp_prob = summarise_imp_prob(I_inc)
    CSV.write(vars.path_imp_fil, df_imp_prob)
end

save_filtered_imp(vars1_IP10)
save_filtered_imp(vars5_IP7)
save_filtered_imp(vars5_IP14)
save_filtered_imp(vars5_IP10)

# ## Generation based Sankey diagram, 0th (Japan), 1st, 2nd, 3rd.

function gen_based_Sankey_diagram(path; nmax=10000)
    path_objs = fetch_sim_paths(path)
    df_sum = DataFrame()

    n = minimum([length(path_objs), nmax])
    @showprogress for i in 1:n
        res = JLD2.load(path_objs[i])["res"]
        global df = DataFrame(res.import_event)
        if size(df)[1] == 0; continue; end

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
                if any(cond); continue; end

                r = filter( x-> x["import"] == df[i, :export], pre_df)[1, :]
                df[i, :gen_index] = r.gen_index + 1
                df[i, :seq_index] = r.seq_index
                pre_df = vcat(pre_df, df[i, :] |> DataFrame)
            end
        end
        df_sum = vcat(df_sum, pre_df)
    end
    return df_sum
end

df_upd = gen_based_Sankey_diagram(vars5_IP10.path_sim_res)
CSV.write(vars5_IP10.path_exp_imp_gen, df_upd)

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

df_fil = filter_exp_gen(vars5_IP10)
CSV.write(vars5_IP10.path_exp_imp_gen_fil, df_fil)

# # Conditional params
# Note: This section may cause an error if individual simulation data are lacking. 

include("model_intercountry.jl")
include("vis_util.jl")

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

    β_med = median(βs)
    β_025 = quantile(βs, 0.025)
    β_975 = quantile(βs, 0.975)

    β_med_r = round(β_med, digits=2)
    β_025_r = round(β_025, digits=2)
    β_975_r = round(β_975, digits=2)
    println("β : $(β_med_r) ($(β_025_r), $(β_975_r))")

    SAR_med = @pipe ( 1- exp(-β_med))*100  |> round(_, digits=3)
    SAR_025 = @pipe ( 1- exp(-β_025))*100  |> round(_, digits=3)
    SAR_975 = @pipe ( 1- exp(-β_975))*100  |> round(_, digits=3)
    println("SAR : $(SAR_med) ($(SAR_025), $(SAR_975))")
end

cut_off = 10
tp = "Korea"

conditional_beta_SAR(vars5_IP10, cut_off, tp; n_sim=50_000)

conditional_beta_SAR(vars1_IP10, cut_off, tp; n_sim=50_000)

conditional_beta_SAR(vars5_IP7, cut_off, tp; n_sim=50_000)

conditional_beta_SAR(vars5_IP14, cut_off, tp; n_sim=50_000)


