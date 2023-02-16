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

include("C:/Workdir/Develop/PF_simulations/My_functions.jl")
include("C:/Workdir/Develop/Hosting_Capacity/HC_functions.jl")
#include("My_ref/My_ref_mn_no_flex.jl")
include("My_ref/My_ref_mn_all_flex.jl")
include("Functions/Mn_functions.jl")
include("Functions/Profiles_functions.jl")
include("Functions/Solve_daily_powerflow.jl")

#=
    Model that performs PF exploiting flexibility in a range of representative days. 
    Demand response and curtailment supported. 
=#

# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//"*file_name
net_data = parse_file(file_path)

######################################
########### PARAMETERS ###############
######################################

parameters = Dict(
    "simulation_periods" => 24,

    "flex" => 50,
    "congestion_limit" =>55,  	             # congestion limit %
    "threshold" => 1,                        # value used to identify branches with current rating higher than threshold
    "curtailment" => 50,                     # Allowed curtailment %

    "multiple_gen" => true,                      # If true, more than one generator is present in the grid
    "size" => 4,                              # size of generator/s
    "seed" => 99,                              # seed for random DG buses choice
    "gen_number_per_feeder" => 2,              # number of random DGs per feeder
    "gen_locations" => ["Ispra", "Ragusa"],    #["Ispra", "Ragusa", "Clamart"], # Generators profiles location
    
    
    "load_locations" => ["Chiara"],             # location of load profiles

)

######################################
########### INPUT DATA ###############
######################################

pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance, n_representative_days, feeder_ID, paths, gen_ID = elaborate_input_data(net_data,parameters)
plot_profiles(pv_profiles,load_profiles,pv_profile_key,load_profile_key, frequency_of_occurance)

AGG_RES = Dict{Int64, Any}()        #Aggregate Results
AGG_NET_DATA = Dict{Int64, Any}()     #Aggregate net data
quick_summary = Dict{Int64, OrderedDict}( r_day => OrderedDict() for r_day in keys(frequency_of_occurance))   #Store main takeaways

for repr_day in keys(frequency_of_occurance) #solve the optimization model for each repr_day, indipendently 
    
    # Prepare model, add gen and load profiles
    mn_net_data, production_hours = prepare_model(
                                                    net_data, 
                                                    repr_day, 
                                                    pv_profiles, 
                                                    pv_profile_key, 
                                                    load_profiles, 
                                                    load_profile_key,
                                                    parameters["simulation_periods"], 
                                                    parameters["gen_locations"],
                                                    parameters["load_locations"], 
                                                    parameters["seed"],
    )

    # Solve model 
    mn_net_data, result = solve_model(mn_net_data, build_mn_pf_DR_DGC) 
    
    # Get informations from results
    flexible_nodes, tot_load, tot_DR, DG_curtailment,DG_production, max_branch_loading, voltage_profile, abs_max_min_volt = evaluate_results(
                                                                                                            mn_net_data, 
                                                                                                            result, 
                                                                                                            production_hours, 
                                                                                                            parameters["threshold"], 
                                                                                                            feeder_ID, 
                                                                                                            paths
    )
    push!(AGG_NET_DATA, repr_day => mn_net_data)
    push!(AGG_RES, repr_day => result)

    #hourly_printing_statements(result, mn_net_data, flexible_nodes, tot_load,DG_curtailment,max_branch_loading, abs_max_min_volt)    
    mn_plot_cumulative_daily_profile(repr_day, parameters["simulation_periods"], mn_net_data, ["nominal generation", "load nominal", "actual generation", "final load"])
    mn_plot_voltage_variation(feeder_ID, voltage_profile, repr_day, "errorband")
    #mn_plot_united_feeder_voltage_profile(voltage_profile, feeder_ID)
    
    update_summary(quick_summary[repr_day], result, DG_production, tot_DR, tot_load, abs_max_min_volt, max_branch_loading, frequency_of_occurance[repr_day])
    
end

println("Demand curtailed in a year (%): ", round(1-sum(x["actual_MWh"]*x["frequency"] for (i,x) in quick_summary)/sum(x["nominal_MWh"]*x["frequency"] for (i,x) in quick_summary), digits = 4)*100)
println("Energy produced curtailed in a year (%): ", round(sum(x["curtailment"]*x["frequency"] for (i,x) in quick_summary)/sum(x["nominal_production_MWh"]*x["frequency"] for (i,x) in quick_summary), digits = 4)*100)