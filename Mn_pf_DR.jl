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
include("My_ref/My_ref_mn_only_DR.jl")
include("Functions/Mn_functions.jl")

#=

    Compute Mutlinetwork PF with only Demand Side Flexibility

=#

#Parameters

simulation_periods = 24

flex = 20                              # flexibility offered %
congestion_limit = 75  	             # congestion limit %
threshold = 20                        # value used to identify branches with current rating higher than threshold


size = 1                              # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 1              # number of random DGs per feeder


# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//Networks//"*file_name
net_data = parse_file(file_path)

load_file = "Test_profiles/Load_profiles.xlsx"  #load profile file
gen_file = "Test_profiles/PV_profile_2016.json"  #generation profile file

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
production_hours = find_production_hours(pv_profile)

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
pm = PowerModels.instantiate_model(mn_net_data, ACPPowerModel, build_mn_pf)
result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])

@assert result["termination_status"] == LOCALLY_SOLVED
println(result["termination_status"])
PowerModels.update_data!(mn_net_data, result["solution"])

# Evaluate flexibility offered
load_variation, tot_load = mn_calc_flexibility_offered(mn_net_data, result)
flexible_nodes = mn_calc_flexible_nodes(load_variation)
tot_DR = calc_tot_DR(flexible_nodes)

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
println( "Maximum brnach loading: $b_max at hour $nw")

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
#mn_printing_statements(result, mn_net_data, flexible_nodes, tot_load)

#mn_plot_united_feeder_voltage_profile(voltage_profile, feeder_ID)

#save_path = "C:/Users/u0152683/Desktop/Networks/Experiments/Multinetwork/"
#=
for (nwid,net) in mn_net_data["nw"]

    plot = mn_plot_grid(net, "p_flex", "pg", "loading", nwid; save_fig = true, save_path = save_path)

end
=#

