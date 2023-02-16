
# INPUT DATA
function elaborate_input_data(net_data::Dict, param::Dict; candidate_gen_nodes = [])

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

    if param["random_gen_assignment"]
        place_random_gens(net_data, param,feeder_ID, gen_number_per_feeder, size, seed; curt = curtailment/100)  #here we add gen randomly
    else
        for node in candidate_gen_nodes
            add_single_generator(net_data, size, node; curt = curtailment/100)  #here we specify ex-ante where to place them
        end
    end

    # Add allowed curtailmentent for each gen
    net_data["curtailment"] = curtailment
    update_data!(net_data, add_allowed_curtailment(net_data, curtailment/100))

    return pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance, n_representative_days, feeder_ID, paths 
end

# PREPARE MODEL BY ASSIGNING LOAD/GEN PROFILES IN EACH DAY  
function prepare_model(data::Dict,repr_day::Int64, pv_profiles::Tuple, pv_profile_key::Dict,load_profiles::Tuple, load_profile_key::Dict)
    
    simulation_periods = parameters["simulation_periods"]
    gen_locations =  parameters["gen_locations"]
    load_locations = parameters["load_locations"]
    seed = parameters["seed"]
    
    #println("\n REPRESENTATIVE DAY $repr_day\n")
    #println("1) Adding generation and load profiles")
    assign_gen_pv_profile(data, gen_locations; seed = seed)
    assign_load_profile(data, load_locations; seed = seed)  

    # REPLICATE DATA
    #println("2) Replicating data")
    mn_net_data = PowerModels.replicate(data,simulation_periods)

    production_hours = mn_find_production_hours(pv_profiles, repr_day)  # productions hours only of PV

    attribute_load_profile_values(mn_net_data, load_profiles, load_profile_key, repr_day)
    attribute_gen_profile_values(mn_net_data, pv_profiles, pv_profile_key, repr_day)

    return mn_net_data, production_hours
end

#SOLVE MODEL
function solve_model(mn_net_data::Dict, model;flag::Bool=false)

    # Solve PF 
    #println("3) SOLVING MODEL")
    pm = PowerModels.instantiate_model(mn_net_data, ACPPowerModel, model) #build_mn_pf_DR_DGC
    result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])

    try
        @assert result["termination_status"] == LOCALLY_SOLVED
        #println(result["termination_status"])
        PowerModels.update_data!(mn_net_data, result["solution"])
        #println("OF :", result["objective"])
        return mn_net_data, result, flag
    catch ex  # if solution is unfeasible, stop iteration
        if ex isa AssertionError 
            println("Solution not feasible")
            flag = true
            return mn_net_data, result, flag
        end
    end

end



# RESULTS INTERPRETATION
function evaluate_results(mn_net_data::Dict, result::Dict, production_hours::Vector{Int64}, param::Dict, feeder_ID::Dict, paths::Vector{Any})

    load_variation, tot_load, flexible_nodes, tot_DR = very_ugly_thing_to_do(4)
    DG_curtailment, DG_production= very_ugly_thing_to_do(2)

    if param["DR"]
        # Evaluate flexibility offered
        load_variation, tot_load = mn_calc_flexibility_offered(mn_net_data, result)
        flexible_nodes = mn_calc_flexible_nodes(load_variation)
        tot_DR = calc_tot_DR(flexible_nodes)
    end

    if param["DGC"]
        # Evaluate curtailment performed
        DG_curtailment = mn_calc_curtailment(mn_net_data, result, production_hours)
        DG_production = mn_calc_tot_production(mn_net_data, result)
    end

    mn_net_data["per_unit"] = true

    # Compute branch flows and losses
    flows = mn_calc_branch_flow(mn_net_data)
    update_data!(mn_net_data, flows)

    attribute_pf_direction(mn_net_data)

    network_losses, tot_network_losses = mn_calc_power_losses(mn_net_data)
    update_data!(mn_net_data, network_losses)

    # Update mn_net_data with branch loadings and provide dict with branch loadings
    branch_loading = mn_calc_branch_loading(mn_net_data)
    max_branch_loading = find_max_branch_loading(branch_loading)

    # Get voltage profiles 
    voltage_profile = calc_voltage_profile(mn_net_data, result, feeder_ID, paths)
    max_min_volt = find_max_min_voltage(voltage_profile)  #max and min voltage for each feeder for each hour
    abs_max_min_volt = find_absolute_max_min_voltage(max_min_volt) # Find the maximum and minimu voltage for each feeder happening in the whole simulation period
    
    return flexible_nodes, tot_load, tot_DR, DG_curtailment, DG_production, max_branch_loading, voltage_profile, abs_max_min_volt
end

function update_summary(summary::OrderedDict, param::Dict, result::Dict, DG_data::Dict, tot_DR::Dict, tot_load, abs_volt::Dict, max_branch_loading::OrderedDict, frequency::Float64,  prod_hours::Vector{Int64})

    summary["OF"] = result["objective"]

    if param["DGC"]
        summary["curtailment"] = sum(x["p_nominal"]- x["p_actual"] for (h,x) in DG_data)  #total curtailment
    end
    summary["actual_production_MWh"] = sum(x["p_actual"] for (h,x) in DG_data)
    summary["nominal_production_MWh"] = sum(x["p_nominal"] for (h,x) in DG_data)
    summary["production_hours"] = prod_hours
    
    if param["DR"]
        summary["up_flex"] = sum(flex["p_flex_up"] for (hour,flex) in tot_DR)
        summary["down_flex"] = sum(flex["p_flex_dwn"] for (hour,flex) in tot_DR)
        summary["nominal_consumption_MWh"] = sum(x["p_demand_nominal"] for (hour,x) in tot_load)
        summary["actual_consumption_MWh"] = sum(x["p_demand_actual"] for (hour,x) in tot_load)
    end

    b_max, nw = findmax(collect(values(max_branch_loading)))
    nw = collect(keys(max_branch_loading))[nw]
    summary["max_branch_loading"] = b_max
    summary["max_branch_loading_hour"] = nw

    d = abs_volt
    vmin = d[reduce((x, y) -> d[x]["vmin"] ≤ d[y]["vmin"] ? x : y, keys(d))]["vmin"]
    vmax = d[reduce((x, y) -> d[x]["vmax"] >= d[y]["vmax"] ? x : y, keys(d))]["vmax"]

    f_ids_max = findall(x -> x["vmax"] == vmax, d)  # all feeders experiencing same max voltage value 
    hrs_max = unique!(reduce(vcat, [d[f_id]["hour_vmax"] for f_id in f_ids_max]))

    f_ids_min = findall(x -> x["vmin"] == vmin, d)  # all feeders experiencing same max voltage value 
    hrs_min = unique!(reduce(vcat, [d[f_id]["hour_vmin"] for f_id in f_ids_min]))

    summary["vmin"] = vmin
    summary["hour_vmin"] = hrs_min
    summary["vmax"] = vmax
    summary["hour_vmax"] = hrs_max

    summary["frequency"] = round(Int64,frequency)


end

function heat_map(net_data::Dict, quick_summary::Dict, AGG_NET_DATA::Dict{Int64,Any}, parameters::Dict; show_plot = false)

    br_year_load = mn_calc_branch_yearly_loading(quick_summary,AGG_NET_DATA)
    yearly_prod_hours = Int.(sum(length(info["production_hours"])*info["frequency"] for (day, info) in quick_summary))
    #mn_plot_branch_yearly_loading(quick_summary, br_year_load,"83")
    branch_rank, branch_loadings = rank_branches(quick_summary, br_year_load, net_data, parameters)

    if show_plot
        plot_grid(net_data,"basic", "basic", "rank";display_flow=false)
    end
    #=
    for branch_id in branch_rank["red"]
        println(branch_id, " => ",sens[6881][branch_id]*100)
    end
    =#
    println( "TOT OF ", sum(x["OF"]*x["frequency"] for (i,x) in quick_summary))
    println()
    for branch_id in branch_rank["red"]
        println("- $branch_id: ", branch_loadings[branch_id]["red_RPF_hours"], " hours")
    end

    return branch_rank, branch_loadings, br_year_load
end