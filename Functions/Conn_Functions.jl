function sort_dict(d::Dict)
    # Sort the keys based on the numbers they have
    sorted_keys = sort(collect(keys(d)), by = x -> parse(Int, filter(isnumeric, x)))
    # Create a new dictionary with the sorted keys
    return  OrderedDict(sorted_keys[i] => d[sorted_keys[i]] for i in 1:length(sorted_keys))
end

# Used to place generators randomly in the grid
function place_random_gens(net_data::Dict, param::Dict, feeder_ID::Dict, gen_number_per_feeder::Int64, size::Int64, seed::Int64; curt::Float64=0.0)

    if param["multiple_gen"]
        #Random choice of buses
        random_generators = get_random_generators(feeder_ID, gen_number_per_feeder, seed)  

        # Add new generators
        add_generators(net_data, random_generators, size; curt = curt)
    else
        
        # Get random generator
        active_node = get_random_DG(seed,net_data)
        
        # Add a random generator
        add_single_generator(net_data, size, active_node; curt = curt)
    end

end

# Add generator with curtailment capabilities 
function add_single_generator(net_data::Dict, size, gen_bus; curt::Float64=0.0)

    if isa(gen_bus, String)
        gen_bus = parse(Int64, gen_bus)
    end
    
    i = length(net_data["gen"]) + 1

    net_data["gen"]["$i"] = Dict("p_nominal" => size, 
    "q_nominal"=>0, 
    "pg" =>size, 
    "qg" =>0, 
    "pmin" => size * (1-curt), 
    "pmax"=>size, 
    "qmin" =>0, 
    "qmax"=>0, 
    "gen_bus" => gen_bus, 
    "gen_status"=>1, 
    "index" => i, 
    "source_id" => ["gen", i],
    "allowed_curt" => curt)
    net_data["bus"]["$gen_bus"]["bus_type"] = 2
    
end

function remove_single_generator(data::Dict, DG_data)

    gen_bus = DG_data["Connection_node"]
    DG = findall(x-> x["gen_bus"] == gen_bus, net_data["gen"])[1]
    delete!(data["gen"], DG)
    
    return data
end

#=
Each generator will get installed a power equal to size_std. 
These generators are added to the model. 
The list of generators are passed as a Dict, where you give a list of buses for each feeders. This is done by calling before get_random_generators
=#
function add_generators(net_data::Dict, generators::Dict{Any, Any}, size_std; curt::Float64=0.0)
    
    x = []
    [append!(x,v) for v in values(generators)]

    for (i,gen) in enumerate(x)
        i+=1
        net_data["gen"]["$i"] = Dict("p_nominal" => size_std, "q_nominal"=>0, "pg" =>size_std, "qg" =>0, "pmin" => size_std * (1-curt) , "pmax"=> size_std , "qmin" =>0, "qmax"=>0, "gen_bus" => gen, "gen_status"=>1, "index" => i, "source_id" => ["gen", i])
        net_data["bus"]["$gen"]["bus_type"] = 2
    end
end

function very_ugly_thing_to_do(n::Int64)
    return [Dict{Any,Any}() for i in 1:n]
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
function mn_plot_branch_yearly_loading(quick_summary::Dict, br_year_load::Dict{Int64,Dict},branch_id::String)

    
    
    d2 = obtain_branch_yearly_loading(quick_summary, br_year_load, branch_id)
     
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
function obtain_branch_yearly_loading(quick_summary::Dict,br_year_load::Dict{Int64,Dict}, branch_id::String)

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
function rank_branches(quick_summary::Dict, br_year_load::Dict, data::Dict, param::Dict)

    branch_loadings = Dict(br_id => Dict() for br_id in keys(data["branch"]))
    branch_rank = Dict{String, Vector{String}}()

    for branch_id in keys(branch_loadings)

        br_loadings = obtain_branch_yearly_loading(quick_summary,br_year_load, branch_id)
        NPF = br_loadings[br_loadings[!,:type].=="NORMAL",:]
        critical_hours_NPF = size(NPF[NPF[!,:loading].>=param["α"],:])[1]
        if critical_hours_NPF>param["β"]
            branch_loadings[branch_id]["red_NPF_hours"] = critical_hours_NPF
            branch_loadings[branch_id]["orange_NPF_hours"] = 0
        else
            branch_loadings[branch_id]["red_NPF_hours"] = 0
            branch_loadings[branch_id]["orange_NPF_hours"] = critical_hours_NPF
        end
        
        RPF = br_loadings[br_loadings[!,:type].=="REVERSE",:]
        critical_hours_RPF = size(RPF[RPF[!,:loading].>=param["α"],:])[1]
        if critical_hours_RPF>param["β"]
            branch_loadings[branch_id]["red_RPF_hours"] = critical_hours_RPF
            branch_loadings[branch_id]["orange_RPF_hours"] = 0
        else
            branch_loadings[branch_id]["red_RPF_hours"] = 0
            branch_loadings[branch_id]["orange_RPF_hours"] = critical_hours_RPF
        end
        #=
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
        =#
    end
    #=
    red_hour_lim = parameters["red_hour_lim"]
    orange_hour_lim = parameters["orange_hour_lim"]

    red_branches = findall(x->sum(x["red_r_hours"]+x["red_n_hours"]) > red_hour_lim,branch_loadings)
    orange_branches = findall(x->sum(x["orange_r_hours"]+x["orange_n_hours"]) > orange_hour_lim,branch_loadings)
    orange_branches = orange_branches[.!(in(red_branches).(orange_branches))]
    non_green_branches = vcat(red_branches,orange_branches)
    =#

    red_branches = findall(x->sum(x["red_RPF_hours"]+x["red_NPF_hours"]) > 0, branch_loadings)
    orange_branches = findall(x->sum(x["orange_RPF_hours"]+x["orange_NPF_hours"]) > 0,branch_loadings)
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

# Provides info about direction of PF in each branch 
function attribute_pf_direction(data::Dict)

    for (nwid, net) in data["nw"]
        [data["nw"][nwid]["branch"][br_id]["flow_dir"] = "reverse" for (br_id,br) in net["branch"] if br["pt"]<0.0]
        [data["nw"][nwid]["branch"][br_id]["flow_dir"] = "normal" for (br_id,br) in net["branch"] if br["pt"]>=0.0]
    end

end

# Add DG reuest into net_data
function introduce_DG_requests(data::Dict, DG::String, DG_data::Dict, param::Dict)
    println("############### Evaluating request for $DG ###################\n")
    println("Addition of $DG in the grid...")
    size = DG_data["P_peak"]
    node = DG_data["Connection_node"]
    curt = param["curtailment"]

    println( "P = ",size, "MWp")
    println( "Node = ",node)

    add_single_generator(data, size, node; curt = curt/100) 
    println("Done!")
end

function update_gen_curt(data::Dict, gen_bus::Int64, curt::Int64)
    gen_id = findall(x -> x["gen_bus"] == gen_bus, data["gen"])[1]
    data["gen"]["$gen_id"]["allowed_curt"] = curt/100
    println()
    println("----- DG curtailment set at: $curt %")
    println()

end