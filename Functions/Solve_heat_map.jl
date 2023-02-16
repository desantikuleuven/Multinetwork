function place_random_gens(net_data::Dict, param::Dict, feeder_ID::Dict, gen_number_per_feeder::Int64, size::Int64, seed::Int64)

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

end


# INPUT DATA
function elaborate_input_data(net_data::Dict, param::Dict; candidate_gen_nodes = [])

    gen_locations = param["gen_locations"]
    load_locations = param["load_locations"]
    congestion_limit = param["congestion_limit"]
    seed = param["seed"]
    gen_number_per_feeder = param["gen_number_per_feeder"]
    size = param["size"]

    # Obtain profiles
    pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance = get_profiles(gen_locations, load_locations)
    n_representative_days = length(frequency_of_occurance)

    # Add congestion capacity for each line
    cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
    update_data!(net_data, cong_cap)

    # Get feeder info 
    feeder_ID, mv_busbar, paths = get_feeder_info(net_data)

    if param["random_gen_assignment"]
        place_random_gens(net_data, param,feeder_ID, gen_number_per_feeder, size, seed)  #here we add gen randomly
    else
        for node in candidate_gen_nodes
            add_single_generator(net_data, size, node)  #here we specify ex-ante where to place them
        end
    end

    # Create dict for DGs (CALL ONLY BEFORE PF)
    gen_ID = get_gen_info(net_data, feeder_ID)

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
function evaluate_results(mn_net_data::Dict, result::Dict, threshold::Int64, feeder_ID::Dict, paths::Vector{Any})

    DG_production = mn_calc_tot_production(mn_net_data, result)

    # Compute branch flows and losses
    flows = mn_calc_branch_flow(mn_net_data)
    update_data!(mn_net_data, flows)

    attribute_pf_direction(mn_net_data)

    network_losses, tot_network_losses = mn_calc_power_losses(mn_net_data)
    update_data!(mn_net_data, network_losses)

    # Update mn_net_data with branch loadings and provide dict with branch loadings
    branch_loading = mn_calc_branch_loading(mn_net_data, threshold)
    max_branch_loading = find_max_branch_loading(branch_loading)

    # Get voltage profiles 
    voltage_profile = calc_voltage_profile(mn_net_data, result, feeder_ID, paths)
    max_min_volt = find_max_min_voltage(voltage_profile)  #max and min voltage for each feeder for each hour
    abs_max_min_volt = find_absolute_max_min_voltage(max_min_volt) # Find the maximum and minimu voltage for each feeder happening in the whole simulation period
    
    return DG_production, max_branch_loading, voltage_profile, abs_max_min_volt
end


function get_summary(summary::OrderedDict, result::Dict, DG_data::Dict, abs_volt::Dict, max_branch_loading::OrderedDict, frequency::Float64, prod_hours::Vector{Int64})

    summary["OF"] = result["objective"]

    summary["actual_production_MWh"] = sum(x["p_actual"] for (h,x) in DG_data)
    summary["nominal_production_MWh"] = sum(x["p_nominal"] for (h,x) in DG_data)
    summary["production_hours"] = prod_hours

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

function attribute_branch_colors(data::Dict, parameters::Dict)

    orange_lower_lim = parameters["orange_lower_lim"]
    orange_upper_lim = parameters["orange_upper_lim"]

    br_load = mn_calc_branch_loading(mn_net_data, 1)

    red_br = Dict{String, Any}()
    orange_br = Dict{String, Any}()
    green_br = Dict{String, Any}()

    for (hour, loading) in br_load
        nwid = string(parse(Int64,hour)+1)

        red_br[nwid] = findall(x-> x> orange_upper_lim, loading)
        orange_br[nwid] = findall(x-> orange_lower_lim <= x <= orange_upper_lim, loading)
        green_br[nwid] = findall(x-> x < orange_lower_lim, loading)
    end

    for (nwid, net) in data["nw"]
        
        [net["branch"][branch_id]["color"] = "red" for branch_id in red_br[nwid]]
        [net["branch"][branch_id]["color"] = "orange" for branch_id in orange_br[nwid]]
        [net["branch"][branch_id]["color"] = "green" for branch_id in green_br[nwid]]

    end

    return red_br, orange_br, green_br

end

function attribute_pf_direction(data::Dict)

    for (nwid, net) in data["nw"]
        [data["nw"][nwid]["branch"][br_id]["flow_dir"] = "reverse" for (br_id,br) in net["branch"] if br["pt"]<0.0]
        [data["nw"][nwid]["branch"][br_id]["flow_dir"] = "normal" for (br_id,br) in net["branch"] if br["pt"]>=0.0]
    end

end

# In order to make this function run, the function attribute_pf_direction() must be called before.
# This function gives a dict with the branch loadings for every hour od the representative days 
function mn_calc_branch_yearly_loading(quick_summary::Dict, AGG_NET_DATA::Dict{Int64,Any})
    
    year_branch_load = Dict{Int64, Dict}(r_day => Dict() for r_day in keys(quick_summary))

    for r_day in keys(quick_summary)

        prd_hours = quick_summary[r_day]["production_hours"]
        daily_branch_load_normal = Dict( br_id => Vector{Float64}() for (br_id,br) in AGG_NET_DATA[r_day]["nw"]["1"]["branch"])
        daily_branch_load_reverse = Dict( br_id => Vector{Float64}() for (br_id,br) in AGG_NET_DATA[r_day]["nw"]["1"]["branch"])

        for (nwid, net) in AGG_NET_DATA[r_day]["nw"]
            hour = parse(Int64,nwid)-1
            if  hour in prd_hours
                
                [push!(daily_branch_load_normal[br_id], br["loading"]) for (br_id,br) in net["branch"] if br["flow_dir"] == "normal"]
                [push!(daily_branch_load_reverse[br_id], br["loading"]) for (br_id,br) in net["branch"] if br["flow_dir"] == "reverse"]      
                    
            end
        end

        push!(year_branch_load[r_day], "normal" => daily_branch_load_normal)
        push!(year_branch_load[r_day], "reverse" => daily_branch_load_reverse)
    end
    return year_branch_load
end

#Plot yearly branch loading 
function mn_plot_branch_yearly_loading(br_year_load::Dict{Int64,Dict},branch_id::String)

    
    
    d2 = obtain_branch_yearly_loading(br_year_load, branch_id)
     
    p1 = d2|> @vlplot(title = "Branch $branch_id yearly loading", 
    width = 1000, height = 600,
    mark = {:line, interpolate= "monotone"}, 
    x = "time:q", 
    y = "loading:q", 
    color = {:type, legend = {title ="PF direction"}}
    ) 
    # Lines for 80% and 50% loading 
    c = DataFrame(lim = [80,50], symb = ["R", "O"])
    p2 = c|> @filter(_.symb == "R") |> @vlplot(width = 1000, height = 600, mark = {:rule, color = "red"}, y = {"lim:q",title = nothing}, size = {value = 2})
    p3 = c|> @filter(_.symb == "O") |> @vlplot(width = 1000, height = 600, mark = {:rule, color = "orange"}, y = {"lim:q",title = nothing}, size = {value = 2})

    
    p = p1 + p2 + p3  # non va se plotto p 

    display(p1)

end

# Returns DataFrame with yearly values of loading for specified branch
function obtain_branch_yearly_loading(br_year_load::Dict{Int64,Dict}, branch_id::String)

    d2 = DataFrame(loading = [], type = [], time = [])
    
    rev_load = sort(reduce(vcat,[repeat(data["reverse"][branch_id],quick_summary[day]["frequency"]) for (day,data) in br_year_load]),rev = true)
    norm_load = sort(reduce(vcat,[repeat(data["normal"][branch_id],quick_summary[day]["frequency"]) for (day,data) in br_year_load]),rev = true)
    if !isempty(norm_load)
        append!(d2[!,"loading"], norm_load )
        append!(d2[!,"type"], reduce(vcat,["NORMAL" for i in 1:length(norm_load)]))
        append!(d2[!,"time"], 1:length(norm_load))
    end

    if !isempty(rev_load)
        append!(d2[!,"loading"], rev_load )
        append!(d2[!,"type"], reduce(vcat,["REVERSE" for i in 1:length(rev_load)]))
        append!(d2[!,"time"], 1:length(rev_load))
    end

    n = d2[d2[!,:type].=="NORMAL",:]
    r = d2[d2[!,:type].=="REVERSE",:]
    return d2
end

# Gives me red, orange and green branches. 
function rank_branches(br_year_load::Dict, data::Dict, parameters::Dict)

    branch_loadings = Dict(br_id => Dict() for br_id in keys(data["branch"]))
    branch_rank = Dict{String, Vector{String}}()

    for branch_id in keys(branch_loadings)

        br_loadings = obtain_branch_yearly_loading(br_year_load, branch_id)
        n = br_loadings[br_loadings[!,:type].=="NORMAL",:]
        red_n_hours = size(n[n[!,:loading].>=parameters["orange_upper_lim"],:])[1]
        orange_n_hours = size( n[parameters["orange_lower_lim"] .<= n[!,:loading] .< parameters["orange_upper_lim"],:])[1]
        n_safe_hours = size(n[n[!,:loading].<parameters["orange_lower_lim"],:])[1]

        branch_loadings[branch_id]["red_n_hours"] = red_n_hours
        branch_loadings[branch_id]["orange_n_hours"] = orange_n_hours
        branch_loadings[branch_id]["n_safe_hours"] = n_safe_hours


        r = br_loadings[br_loadings[!,:type].=="REVERSE",:]
        red_r_hours = size(r[r[!,:loading].>=parameters["orange_upper_lim"],:])[1]
        orange_r_hours = size(r[parameters["orange_upper_lim"] .> r[!,:loading].>=parameters["orange_lower_lim"],:])[1]
        r_safe_hours = size(r[r[!,:loading].<parameters["orange_lower_lim"],:])[1]

        branch_loadings[branch_id]["red_r_hours"] = red_r_hours
        branch_loadings[branch_id]["orange_r_hours"] = orange_r_hours
        branch_loadings[branch_id]["r_safe_hours"] = r_safe_hours
    end

    red_hour_lim = parameters["red_hour_lim"]
    orange_hour_lim = parameters["orange_hour_lim"]

    red_branches = findall(x->sum(x["red_r_hours"]+x["red_n_hours"]) > red_hour_lim,branch_loadings)
    orange_branches = findall(x->sum(x["orange_r_hours"]+x["orange_n_hours"]) > orange_hour_lim,branch_loadings)
    orange_branches = orange_branches[.!(in(red_branches).(orange_branches))]
    non_green_branches = vcat(red_branches,orange_branches)

    [data["branch"][i]["rank"] = "red" for i in red_branches]
    [data["branch"][i]["rank"] = "orange" for i in orange_branches]
    [data["branch"][i]["rank"] = "green" for (i,bla) in data["branch"] if !(i in non_green_branches)]
    branch_rank["red"] = red_branches
    branch_rank["orange"] = orange_branches

    return branch_rank, branch_loadings

end

# Compute the branch sensitivity for a variation in the power injection in a specified node 
function calc_branch_ptdf(data::Dict, bus::Vector; var::Float64=-0.2)

    function sensitivity(net::Dict, res_1::Dict, var::Float64)  #sensitivity obtained by adding gen. Default variation in gen is 20%.

        for (gen_id,gen) in net["gen"]
            if gen_id != "1"
                net["gen"][gen_id]["pg"] *= (1+var)
                net["gen"][gen_id]["pmin"] *= (1+var)
                net["gen"][gen_id]["pmax"] *= (1+var)
                net["gen"][gen_id]["p_nominal"] *= (1+var)
            end
        end

        pm = PowerModels.instantiate_model(net, ACPPowerModel, no_flex_with_DG) #build_mn_pf_DR_DGC
        res_2 = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])

        sens = Dict{String, Float64}()

        for (branch_id, branch) in net["branch"]

            # Ref values
            S_ref = branch["rate_a"]

            # initial values
            br = res_1["solution"]["branch"][branch_id]
            S_start = max(abs(complex(br["pt"],br["qt"])) , abs(complex(br["pf"],br["qf"])))

            # final values
            br_2 = res_2["solution"]["branch"][branch_id]
            S_final = max(abs(complex(br_2["pt"],br_2["qt"])) , abs(complex(br_2["pf"],br_2["qf"])))

            #sensitivity 
            ΔS = (S_final - S_start)/S_ref
            ΔP = -var  # P_ref is 1, I prefere to keep denominator as positive

            sens[branch_id] = ΔS/ΔP
        end

        return filter(x->abs(last(x))>0.05,sens)
    
    end

    ptdf = Dict{Int64,Dict}()
    network = deepcopy(data)
    # clean all generators present
    if length(keys(network["gen"]))>1
        [delete!(network["gen"], i) for i in keys(data["gen"]) if i != "1"]
    end

    for bus_id in bus
        
        net = deepcopy(network)
        add_single_generator(net, 1 , bus_id)  # we add 1 MW gen
        
        pm = PowerModels.instantiate_model(net, ACPPowerModel, no_flex_with_DG) #build_mn_pf_DR_DGC
        res_1 = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])
        
        push!(ptdf, bus_id  => sensitivity(net, res_1, var))

    end

    return ptdf
end