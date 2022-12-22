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

include("C:/Workdir/Develop/Networks/PF_simulations/My_functions.jl")
include("C:/Workdir/Develop/Networks/Hosting_Capacity/HC_functions.jl")
include("My_ref/My_ref_mn_no_flex.jl")
include("Functions/Mn_functions.jl")

#=

    Compute Mutlinetwork PF with only Demand Side Flexibility

=#

#Parameters

simulation_periods = 24

congestion_limit = 100  	             # congestion limit %
threshold = 20                        # value used to identify branches with current rating higher than threshold


size = 2                              # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 1              # number of random DGs per feeder


# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//Networks//"*file_name
net_data = parse_file(file_path)

load_file = "Test_profiles/Load_profiles.xlsx"  #load profile file
gen_file = "Test_profiles/PV_profile_2016.json"  #generation profile file

# Add congestion capacity for each line
cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Get feeder info 
feeder_ID, mv_busbar, paths = get_feeder_info(net_data)

#Random choice of buses
random_generators = get_random_generators(feeder_ID, gen_number_per_feeder, seed)  

# Add new generators
add_generators(net_data, random_generators, size)

# Get random generator
active_node = get_random_DG(seed,net_data)

# Add a random generator
#add_single_generator(net_data, size, active_node)

# Create dict for DGs (CALL ONLY BEFORE PF)
gen_ID = get_gen_info(net_data, feeder_ID)

# Choose for each load their profile
choose_load_profiles(net_data, load_file; seed_value = seed)   # ALWAYS CALL BEFORE replicate()

# Replicate data
mn_net_data = PowerModels.replicate(net_data,simulation_periods)

pv_profile = retrieve_gen_profile(gen_file, "summer")
production_hours = find_production_hours(pv_profile)  # productions hours only of PV

#Assign production profile to gens
for (n,net) in mn_net_data["nw"]
    for (i,gen) in net["gen"]
        if i!= 1
            k = parse(Int64,n)
            mn_net_data["nw"][n]["gen"][i]["pg"] *= pv_profile[k]
            mn_net_data["nw"][n]["gen"][i]["qg"] *=  pv_profile[k]
        end
        
    end
end

assign_load_profile(mn_net_data, load_file)
#assign_gen_profile(mn_net_data, load_file, "summer")  #doesn't work, there is a bug. However fnuction does the job correctly

# Solve PF 
pm = PowerModels.instantiate_model(mn_net_data, ACPPowerModel, no_flex_with_DG)
result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])

@assert result["termination_status"] == LOCALLY_SOLVED
println(result["termination_status"])
PowerModels.update_data!(mn_net_data, result["solution"])


# Compute branch flows and losses
flows = mn_calc_branch_flow(mn_net_data)
update_data!(mn_net_data, flows)

network_losses, tot_network_losses = mn_calc_power_losses(mn_net_data)
update_data!(mn_net_data, network_losses)

# Update mn_net_data with branch loadings and provide dict with branch loadings
branch_loading = mn_calc_branch_loading(mn_net_data, 1)
max_branch_loading = find_max_branch_loading(branch_loading)  #gets you the max loading in a day 

# Get voltage profiles 
voltage_profile = calc_voltage_profile(mn_net_data, result, feeder_ID, paths)
max_min_volt = find_max_min_voltage(voltage_profile)  #max and min voltage for each feeder for each hour
abs_max_min_volt = find_absolute_max_min_voltage(max_min_volt) # Find the maximum and minimu voltage for each feeder happening in the whole simulation period

# Identify th enumber of hours each branch is in a red, orange or green status 
branch_zones = mn_calc_branch_zones(net_data, branch_loading)

critical_daylight_hour = critical_hour_max_gen_min_load(mn_net_data)   # hour with max gen, min load
#critical_night_hour = critical_hour_min_gen_max_load(mn_net_data)    # hours with min gen, max load 


r_cong_h, o_cong_h = most_congested_hours(branch_zones)

#=function attribute_load_map(branch_zones::Dict, data::Dict)
    
    for (branch_id, info) in branch_zones



    end
end
=#
most_cong_h = find_most_congested_hour(r_cong_h, o_cong_h, simulation_periods)



#mn_plot_cumulative_daily_profile(simulation_periods, mn_net_data, ["final generation", "load nominal"])