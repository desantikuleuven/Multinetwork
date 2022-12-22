using PowerModels
using PowerModelsAnalytics
using PowerPlots
using JuMP, Ipopt
using Setfield
import InfrastructureModels; const _IM = InfrastructureModels
using DataFrames
using DataStructures
using StatsBase
using Random
using YAML
using Dates
using VegaLite

include("C:/Workdir/Develop/Networks/PF_simulations/My_functions.jl")
include("C:/Workdir/Develop/Networks/Hosting_Capacity/HC_functions.jl")
#include("My_ref/My_ref_mn_no_flex.jl")
include("My_ref/My_ref_mn_all_flex.jl")
include("Functions/Mn_functions.jl")
include("Functions/Profiles_functions.jl")
include("Functions/Solve_daily_powerflow.jl")

#=
    Compute Mutlinetwork PF with only Demand Side Flexibility
=#

# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//Networks//"*file_name
net_data = parse_file(file_path)

######################################
########### PARAMETERS ###############
######################################

parameters = Dict(
    "simulation_periods" => 24,

    "flex" => 50,
    "congestion_limit" =>100,  	             # congestion limit %
    "threshold" => 1,                        # value used to identify branches with current rating higher than threshold
    "curtailment" => 50,                     # Allowed curtailment %

    "multiple_gen" => true,                      # If true, more than one generator is present in the grid
    "size" => 2,                              # size of generator/s
    "seed" => 99,                              # seed for random DG buses choice
    "gen_number_per_feeder" => 2,              # number of random DGs per feeder
    "gen_locations" => ["Ispra", "Ragusa"],#["Ispra", "Ragusa", "Clamart"], # Generators profiles location
    
    
    "load_locations" => ["Chiara"],             # location of load profiles

)

######################################
########### INPUT DATA ###############
######################################

pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance, n_representative_days, feeder_ID, paths, gen_ID = elaborate_input_data(net_data,parameters)
#plot_profiles(pv_profiles,load_profiles,pv_profile_key,load_profile_key, frequency_of_occurance; show_gen = true)

AGG_RES = tuple((DataFrames.DataFrame() for day in 1:n_representative_days)...)  #Aggregate Results
AGG_NET_DATA = tuple((DataFrames.DataFrame() for day in 1:n_representative_days)...)  #Aggregate net data

model = "build_mn_pf_DR_DGC"
#model = "no_flex_with_DG"   # model without any kind of flexibility

quick_summary = Dict{Int64, Any}(i => Dict() for i in keys(frequency_of_occurance))

for (idx,repr_day) in enumerate(keys(frequency_of_occurance)) #solve the optimization model for each repr_day, indipendently 
    
    # Prepare model, add gen and load profiles
    mn_net_data, production_hours = prepare_model(
                                                    net_data, 
                                                    repr_day, 
                                                    pv_profiles, 
                                                    pv_profile_key, 
                                                    load_profiles, 
                                                    load_profile_key,
                                                    simulation_periods, 
                                                    parameters["gen_locations"],
                                                    parameters["load_locations"], 
                                                    seed
    )

    # Solve model 
    mn_net_data, result = solve_model(mn_net_data, build_mn_pf_DR_DGC)
    
    # Get informations from results
    flexible_nodes, tot_load, tot_DR, DG_curtailment, max_branch_loading, voltage_profile, abs_max_min_volt = evaluate_results(
                                                                                                            mn_net_data, 
                                                                                                            result, 
                                                                                                            production_hours, 
                                                                                                            parameters["threshold"], 
                                                                                                            feeder_ID, 
                                                                                                            paths
    )
    append!(AGG_NET_DATA[idx], mn_net_data)
    append!(AGG_RES[idx], result)

    #=
    
    for (nwid, flex) in tot_DR
        if flex["p_flex_up"] != 0
            println("Upwards active flexibility used: ",flex["p_flex_up"],"MW at $nwid")
            quick_summary[repr_day]["Up_flex"] = flex["p_flex_up"]
        elseif flex["p_flex_dwn"] != 0
            println("Dwonwards active flexibility used: ",flex["p_flex_dwn"]," at $nwid")
            quick_summary[repr_day]["Dwon_flex"] = flex["p_flex_dwn"]
        end
    end
    =#

    #hourly_printing_statements(result, mn_net_data, flexible_nodes, tot_load,DG_curtailment,max_branch_loading, abs_max_min_volt)    
    mn_plot_cumulative_daily_profile(repr_day, simulation_periods, mn_net_data, ["nominal generation", "load nominal", "actual generation", "final load"])

    update_summary(quick_summary[repr_day],result,mn_net_data,tot_DR)
    

    
end
