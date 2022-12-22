
# INPUT DATA
function elaborate_input_data(net_data::Dict, param::Dict)

    gen_locations = param["gen_locations"]
    load_locations = param["load_locations"]
    congestion_limit = param["congestion_limit"]
    seed = param["seed"]
    gen_number_per_feeder = param["gen_number_per_feeder"]
    size = param["size"]
    flex = param["flex"]
    curtailment = param["curtailment"]

    # Obtain profiles
    pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance = get_profiles(gen_locations, load_locations)
    n_representative_days = length(frequency_of_occurance)

    # Add flexibility % that each load can offer
    net_data["flex"] = flex
    update_data!(net_data, add_load_flexibility(net_data, net_data["flex"]/100))

    # Add congestion capacity for each line
    cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
    update_data!(net_data, cong_cap)

    # Get feeder info 
    feeder_ID, mv_busbar, paths = get_feeder_info(net_data)

    if param["multiple_gen"]
        #Random choice of buses
        random_generators = get_random_generators(feeder_ID, gen_number_per_feeder, seed)  

        # Add new generators
        add_generators(net_data, random_generators, size)
    else

        # Get random generator
        active_node = get_random_DG(seed,net_data)

        # Add a random generator
        add_single_generator(net_data, size, active_node)
    end

    # Create dict for DGs (CALL ONLY BEFORE PF)
    gen_ID = get_gen_info(net_data, feeder_ID)

    # Add allowed curtailmentent for each gen
    net_data["curtailment"] = curtailment
    update_data!(net_data, add_allowed_curtailment(net_data, curtailment/100))

    return pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance, n_representative_days, feeder_ID, paths, gen_ID 
end

# MODEL 
function prepare_model(data::Dict,repr_day::Int64, pv_profiles::Tuple, pv_profile_key::Dict,load_profiles::Tuple, load_profile_key::Dict, simulation_periods::Int64, gen_locations::Vector{String}, load_locations::Vector{String}, seed::Int64)
    println("\n REPRESENTATIVE DAY $repr_day\n")
    println("1) Adding generation and load profiles")
    assign_gen_pv_profile(data, gen_locations; seed = seed)
    assign_load_profile(data, load_locations; seed = seed)  

    # REPLICATE DATA
    println("2) Replicating data")
    mn_net_data = PowerModels.replicate(data,simulation_periods)

    production_hours = mn_find_production_hours(pv_profiles, repr_day)  # productions hours only of PV

    attribute_load_profile_values(mn_net_data, load_profiles, load_profile_key, repr_day)
    attribute_gen_profile_values(mn_net_data, pv_profiles, pv_profile_key, repr_day)

    return mn_net_data, production_hours
end

function solve_model(mn_net_data::Dict, model)
    # Solve PF 
    println("3) SOLVING MODEL")
    pm = PowerModels.instantiate_model(mn_net_data, ACPPowerModel, model) #build_mn_pf_DR_DGC
    result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])

    @assert result["termination_status"] == LOCALLY_SOLVED
    println(result["termination_status"])
    PowerModels.update_data!(mn_net_data, result["solution"])
    println("OF :", result["objective"])

    return mn_net_data, result

end

# RESULTS INTERPRETATION
function evaluate_results(mn_net_data::Dict, result::Dict, production_hours::Vector{Int64}, threshold::Int64, feeder_ID::Dict, paths::Vector{Any})
    
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
    
    return flexible_nodes, tot_load, tot_DR, DG_curtailment, max_branch_loading, voltage_profile, abs_max_min_volt
end

