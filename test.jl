using PowerModels
using PowerModelsAnalytics
using PowerPlots
using ColorSchemes
using JuMP, Ipopt
using Setfield
using Plots
using VegaLite
import InfrastructureModels; const _IM = InfrastructureModels
using DataFrames
using StatsBase
using Random

include("C:/Workdir/Develop/Networks/PF_simulations/My_functions.jl")
include("C:/Workdir/Develop/Networks/Hosting_Capacity/HC_functions.jl")
include("My_ref/My_ref_mn_all_flex.jl")
include("Functions/Mn_functions.jl")

#=

    Compute Mutlinetwork PF with both DR and DG curt

=#

#Parameters

simulation_periods = 24

flex = 20                              # flexibility offered %
congestion_limit = 75  	             # congestion limit %
threshold = 20                       # value used to identify branches with current rating higher than threshold
curtailment = 20                       # DG curtailment %

size = 9                               # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 1              # number of random DGs per feeder


# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//Networks//"*file_name
net_data = parse_file(file_path)

load_file = "Test_profiles/Load_profiles.xlsx"  #load profile file
gen_file = "Test_profiles//PV_profile_2016.json"  #generation profile file

# Add flexibility % that each load can offer
net_data["flex"] = flex
update_data!(net_data, add_load_flexibility(net_data, net_data["flex"]/100))

# Add congestion capacity for each line
cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Get feeder info 
feeder_ID, mv_busbar, paths = get_feeder_info(net_data)

#Random choice of buses
random_generators = get_random_generators(feeder_ID, gen_number_per_feeder, seed)  

# Add new generators
add_curtailable_generators(net_data, random_generators, size, curtailment/100)

# Get random generator
active_node = get_random_DG(seed,net_data)

# Add a random generator
#add_single_curtailable_generator(net_data, size, active_node; curt = curtailment/100)

# Create dict for DGs (CALL ONLY BEFORE PF)
gen_ID = get_gen_info(net_data, feeder_ID)

# Add allowed curtailmentent for each gen
net_data["curtailment"] = curtailment
update_data!(net_data, add_allowed_curtailment(net_data, curtailment/100))

# Choose for each load their profile
choose_load_profiles(net_data, load_file; seed_value = seed)   # ALWAYS CALL BEFORE replicate()

# Replicate data
mn_net_data = PowerModels.replicate(net_data,simulation_periods)

# Assign production and consumption profiles 
pv_profile = retrieve_gen_profile(gen_file, "winter")
production_hours = find_production_hours(pv_profile)

for (n,net) in mn_net_data["nw"] #Assign production profile to gens
    for (i,gen) in net["gen"]
        if i!= "1"
            k = parse(Int64,n)
            mn_net_data["nw"][n]["gen"][i]["pg"] *= pv_profile[k]
            mn_net_data["nw"][n]["gen"][i]["qg"] *= pv_profile[k]

            mn_net_data["nw"][n]["gen"][i]["pg_ref"] = mn_net_data["nw"][n]["gen"][i]["pg"]  # used as reference value for curtailment evaluation
            mn_net_data["nw"][n]["gen"][i]["qg_ref"] = mn_net_data["nw"][n]["gen"][i]["qg"]
        end
        
    end
end
assign_load_profile(mn_net_data, load_file)

#assign_gen_profile(mn_net_data, load_file, "summer")  #doesn't work, there is a bug. However fnuction does the job correctly

# Solve PF 
pm = PowerModels.instantiate_model(mn_net_data, ACPPowerModel, build_mn_pf_DR_DGC)
result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])

@assert result["termination_status"] == LOCALLY_SOLVED
println(result["termination_status"])
PowerModels.update_data!(mn_net_data, result["solution"])
println("OF :", result["objective"])

# Evaluate flexibility offered
load_variation, tot_load = mn_calc_flexibility_offered(mn_net_data, result)
flexible_nodes = mn_calc_flexible_nodes(load_variation)
tot_DR = calc_tot_DR(flexible_nodes)

# Evaluate curtailment performed
DG_curtailment = mn_calc_curtailment(mn_net_data, result, production_hours)
DG_production = mn_calc_tot_production(mn_net_data, result)

# Compute branch flows and losses
flows = mn_calc_branch_flow(mn_net_data)
update_data!(mn_net_data, flows)

network_losses, tot_network_losses = mn_calc_power_losses(mn_net_data)
update_data!(mn_net_data, network_losses)

# Update mn_net_data with branch loadings and provide dict with branch loadings
branch_loading = mn_calc_branch_loading(mn_net_data, threshold)
max_branch_loading = find_max_branch_loading(branch_loading)

# Get voltage profiles 
voltage_profile = calc_voltage_profile(mn_net_data, result, feeder_ID, paths)
max_min_volt = find_max_min_voltage(voltage_profile)  #max and min voltage for each feeder for each hour
abs_max_min_volt = find_absolute_max_min_voltage(max_min_volt) # Find the maximum and minimu voltage for each feeder happening in the whole simulation period

b_max, nw = findmax(collect(values(max_branch_loading)))
nw = collect(keys(max_branch_loading))[nw]
println( "Maximum branch loading: $b_max at hour $nw")

for (f_id, values) in abs_max_min_volt
    println()
    println("Feeder $f_id: V max = ",values["vmax"]," at ",values["hour_vmax"],"  V min = ",values["vmin"]," at ",values["hour_vmin"])

end

for (nwid, flex) in tot_DR
    if flex["p_flex_up"] != 0
        println("Upwards active flexibility used: ",flex["p_flex_up"],"MW at $nwid")
    elseif flex["p_flex_dwn"] != 0
        println("Dwonwards active flexibility used: ",flex["p_flex_dwn"]," at $nwid")
    end
end

#mn_printing_statements(result, mn_net_data, flexible_nodes, tot_load,DG_curtailment)


#mn_plot_united_feeder_voltage_profile(voltage_profile, feeder_ID)
#=
for (nwid,net) in mn_net_data["nw"]

    plot = mn_plot_grid(net, "p_flex", "pg", "loading", nwid;save_fig = true, save_path = "C:/Users/u0152683/Desktop/Networks/Experiments/Multinetwork/")

end
=#
#plot = mn_plot_grid(mn_net_data["nw"]["14"], "p_flex", "pg", "loading", "14")
function aliases( value::String)

    value = lowercase(value)

    if value in ["load final", "demand final", "load actual", "demand actual","load_p", "final load", "final demand"]
        return "load_actual"

    elseif value in ["load nominal", "demand nominal", "load", "demand", "pd", "nominal load", "nominal demand"]
        return "load_nominal"

    elseif value in ["final generation", "final production", "actual generation", "actual production","pg actual", "pg final"]
        return "pg_actual"

    elseif value in ["nominal generation", "nominal production", "generation", "production", "pg", "pg nominal", "pg ref"]
        return "pg_nominal"
    else
        println("Proprieties to display invalid")
    end
end

function calc_hourly_stuff(df::PowerModelsDataFrame, value::String)

    if value == "load_actual"
        return sum(df.load.load_p)

    elseif value == "load_nominal"
        return sum(df.load.pd)

    elseif value == "pg_actual"
        return sum(df.gen.pg[df.gen.gen_bus.!=1])

    elseif value == "pg_nominal"
        return sum(df.gen.pg_ref[df.gen.gen_bus.!=1])
    end
end

function update_columns(df::DataFrame, PMD::PowerModelsDataFrame, prop_to_show::Vector{String}, idx::Int64)

    for prop in prop_to_show
        df[idx, prop] = calc_hourly_stuff(PMD, prop)
    end

end

function convert_to_single_df(df::DataFrame, sim_periods::Int)
    power = Vector{Float64}()
    type = []
    for col_name in filter(x -> x!= "hour", names(df))
        power = vcat(power, df[!,col_name])
        type =  vcat(type,reduce(vcat, [col_name for i in 1:sim_periods]))
    end

    hours = repeat(1:sim_periods,length(filter(x -> x!= "hour", names(df))))

    return DataFrame(power = power, type = type, hour = hours)
end

function mn_plot_cumulative_daily_profile(sim_periods::Int64, data::Dict, prop_to_show::Vector{String})

    df = DataFrame(hour = 1:sim_periods)  #create dataframe where values for each hour will be stored
    proprieties_to_show = Vector{String}()

    for prop in prop_to_show
        prop = aliases(prop)
        println(prop)
        push!(proprieties_to_show, prop)
        df[!, prop] = zeros(sim_periods)
    end

    for (nwid, net) in data["nw"]
    
        PMD = PowerModelsDataFrame(net)
        idx = parse(Int64,nwid)

        update_columns(df, PMD, proprieties_to_show, idx)
    end
    
    data_df = convert_to_single_df(df, sim_periods)
    
    data_df |> 
    @vlplot(
        :line, 
        color = {
            :type,
            legend = { title = "Line type"}
        },
        width = 500, height = 300,
        x = {
            :hour,
            title = "Hours"
        },

        y = {
            :power,
            title = "Power (MW)"
        },
    
    )  
    
end

#mn_plot_cumulative_daily_profile(simulation_periods, mn_net_data, ["nominal generation", "final generation", "nominal load", "final load"])
mn_plot_cumulative_daily_profile(simulation_periods, mn_net_data, ["nominal generation", "nominal load"])