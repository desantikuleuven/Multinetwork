using DataFrames
using XLSX
using JSON
using DataStructures
using Dates
using PowerPlots
using PowerPlots.Experimental
using ColorSchemes
using Plots

# Assign the load profile to each load in the multinetwork network 
function assign_load_profile(data::Dict, load_file::String)

    @assert data["multinetwork"] == true

    load_profiles = retrieve_load_profiles(load_file)

    #@assert length(data["multinetwork"]) == size()

    for (n,net) in data["nw"]
        for (i,load) in net["load"]
            k = parse(Int64,n)
            data["nw"][n]["load"][i]["pd"] = load["pd"] * load_profiles[load["load_profile"]][k]
            data["nw"][n]["load"][i]["qd"] = load["qd"] * load_profiles[load["load_profile"]][k]

            data["nw"][n]["load"][i]["pd_ref"] = data["nw"][n]["load"][i]["pd"]
            data["nw"][n]["load"][i]["qd_ref"] = data["nw"][n]["load"][i]["qd"]
        end
    end


end

function choose_gen_pv_profile(data::Dict, profiles::DataFrame; seed_value::Int)
    Random.seed!(seed_value)
    profile_names = names(profiles)
    deleteat!(profile_names, findall(x->x=="time",profile_names))
    
    for (gen_id, gen) in data["gen"]
        if gen_id!= "1"
            data["gen"][gen_id]["gen_profile"] = sample(profile_names,1)[1]  #random assingment of profile 
        end
    end

end

# Attribute to each load the profil desired, mode random
function choose_load_profiles(data::Dict, load_file::String; seed_value::Int)

    load_profiles = retrieve_load_profiles(load_file)

    Random.seed!(seed_value)
    
    for (i,load) in data["load"]
        data["load"][i]["load_profile"] = sample(collect(keys(load_profiles)),1)[1]
    end

end

# Same function as calc_voltage_profile but without showing voltage profiles, needed just to get feeders info
function get_feeder_data(net_data::Dict, file_name)
    
    # INPUT FOR PATH CALC

    down_node, g, map = downstreamcalcs(net_data)

    extremes = []
    [push!(extremes, i) for (i,j) in down_node if j ==[]]

    tot_paths = [Tuple[]]
    dist = AbstractFloat[]

    len = Dict((map[j["from_b"]],map[j["to_b"]]) => j["length"] for (ind,j) in net_data["distance"])

    inv_map = inverse_map(net_data)
    node_dist = []

    tot_paths, dist, node_dist = paths_calc(g, extremes, len, dist, tot_paths, inv_map, file_name, node_dist)

    paths = []
    feeders = []

    mv_busbar = 0

    #COUNT NUMBEER OF FEEDERS
    for (path,line) in enumerate(tot_paths)
    
        if line[1][2]!=1 && line[1][2] ∉ feeders
            push!(feeders, line[1][2])
        end
    end

    # CREATE DICT CONTAINING ALL DATA RELEVANT FOR FEEDERS
    alphabet = 'A':'Z'
    feeder_ID = Dict( j => Dict("Name" => "Feeder $i", "Paths" => [], "Paths_ID" => [], "Paths_distance"=> []) for (i,j) in zip('A':alphabet[length(feeders)], feeders) )

    # Write tot_paths in a vectorial way. 
    for vec in tot_paths

        feed = []
        for j in 1:length(vec)

            push!(feed,vec[j][1])
            
        end
    
        push!(paths,feed)  #Rewrite how paths are made

    end

    for (idx,path) in enumerate(paths)

        if path[end]!=1 && path[2] in keys(feeder_ID)
            push!(feeder_ID[path[2]]["Paths"], path)
            push!(feeder_ID[path[2]]["Paths_ID"], idx)
            push!(feeder_ID[path[2]]["Paths_distance"], node_dist[idx])
        end

    end
    
    # Add ID of branches belonging to the same feeder

    for (ref, data) in feeder_ID

        branches = []
    
        for path in feeder_ID[ref]["Paths"]
            
            for (ind,j) in enumerate(path)
        
                if ind<length(path)
        
                    t_bus = path[ind]
                    f_bus = path[ind+1]
        
                    for (idx,branch) in net_data["branch"]
        
                        if f_bus == branch["f_bus"] && t_bus == branch["t_bus"]
                            if idx ∉ branches
                                push!(branches,idx)
                            end
                        end
                    end
                end
            end
        end
    
        feeder_ID[ref]["Branches"] = branches
    
    end
    
    # Add ID of bus belonging to the same feeder
    for (ref, data) in feeder_ID

        nodes = []

        for path in feeder_ID[ref]["Paths"]
            [push!(nodes, bus) for bus in path if bus ∉ nodes]
        end

        feeder_ID[ref]["Buses"] = nodes
        mv_busbar = nodes[1]
    end

    return feeder_ID, mv_busbar
end

# Assign the generation profile to each generator in the multinetwork network 
function assign_gen_profile(data::Dict, gen_file::String, season::String)

    @assert data["multinetwork"] == true

    pv_profile = retrieve_gen_profile(gen_file, season)

    for (n,net) in data["nw"]
        for (i,gen) in net["gen"]
            if i!= "1"
                k = parse(Int64,n)
                data["nw"][n]["gen"][i]["p_nominal"] *= pv_profile[k]
                data["nw"][n]["gen"][i]["q_nominal"] *=  pv_profile[k]
            end
            
        end
    end


end

# Get the load profiles from the Excel sheet
function retrieve_load_profiles(load_file)
    xf = XLSX.readxlsx(load_file)
    profile_types = XLSX.sheetnames(xf)

    profiles = Dict{String,Dict{Int64,Float64}}(profile => Dict() for profile in profile_types) 

    for sheet in profile_types

        data = DataFrame(XLSX.readtable(load_file,sheet))
    
        for (i,val) in enumerate(data[:,2])  #second column of data is where values are stored
            profiles[sheet][i] = val
        end

    end

    return profiles
end

# From the yearly generation profile, obtain the monthly profile specified in the argument. 
function monthly_gen_profile(pv_prod::OrderedDict, ref_year::Int, month::Int)

    # Identify relevent month hours
    if month==12
        period = collect(Dates.DateTime(ref_year, month):Dates.Hour(1):Dates.DateTime(ref_year+1, 1))
        pop!(period)
    else
        period = collect(Dates.DateTime(ref_year, month):Dates.Hour(1):Dates.DateTime(ref_year, month+1)) # period contains all the hours in the ref year and month
        pop!(period)
    end

    # Obtain monthly profile of pv, separate it for each day
    monthly_profile = OrderedDict{Int64, Vector{Float64}}(i => Float64[] for i in sort(Dates.day.(period)))

    for date in period

        day = Dates.day(date) # day of the month 
        push!(monthly_profile[day], pv_prod[date])

    end

    return monthly_profile

end

# gives back PV generation profile during summer or winter.
# Summer profile is obtained by giving the day corresponding to highest peak power during the whole summer
# Winter profile is obtained by giving the day corresponding to lowest peak power during the whole winter
function retrieve_gen_profile(gen_file::String, season::String="summer")
    data = JSON.parsefile(gen_file)

    ref_year = data["inputs"]["meteo_data"]["year_max"]
    p_ref = data["inputs"]["pv_module"]["peak_power"]  #kWp

    ref_dates = collect(Dates.DateTime(ref_year):Dates.Hour(1):Dates.DateTime(ref_year+1))
    pop!(ref_dates)

    @assert length(ref_dates) == length(data["outputs"]["hourly"])

    pv_prod = OrderedDict{DateTime, Float64}(date => (pv["P"]/p_ref)*10^-3 for (date,pv) in zip(ref_dates, data["outputs"]["hourly"]))

    if season in ["summer", "Summer", "estate", "Estate", "optimistic", "Optimistic"]

        ref = Dict{String, Any}("peak_power" =>0.0 , "mean_daily_power_production" =>0.0 )
        
        for i in [4,5,6,7,8,9,10]  #Evaluate production during summer 

            monthly_profile = monthly_gen_profile(pv_prod, ref_year, i)

            summer_criteria(monthly_profile, i, ref)

        end
        
        println("Choosen summer profile sorresponding to day ", ref["day"],"/",ref["month"], " with peak power of ", ref["peak_power"])
        println("and highest daily average power of ", ref["mean_daily_power_production"], "\n")
     

        daily_profile = monthly_gen_profile(pv_prod, ref_year, ref["month"])[ref["day"]]

        return daily_profile

    elseif season in ["winter", "Winter", "inverno", "Inverno", "pessimistic", "Pessimistic"]

        ref = Dict{String, Any}("peak_power" =>1.0, "mean_daily_power_production" =>1.0)

        for i in [10,11,12,1,2,3]  #Evaluate production during winter 

            monthly_profile = monthly_gen_profile(pv_prod, ref_year, i)

           

            winter_criteria(monthly_profile, i, ref)

        end

        println("Choosen winter profile sorresponding to day ", ref["day"],"/",ref["month"], " with peak power of ", ref["peak_power"],"\n")
        println("and lowest daily average power of ", ref["mean_daily_power_production"])
     

        daily_profile = monthly_gen_profile(pv_prod, ref_year, ref["month"])[ref["day"]]

        return daily_profile
    end

    
end

# Criteria for the selection of the pv generation profile during summer (i.e. most productive)
# Default is based on highest daily average production, other optional criteria is maximum peak power
function summer_criteria(monthly_profile::OrderedDict, month::Int64, ref::Dict)

    daily_avg_prod, daily_p_max = calc_daily_data(monthly_profile)

    highest_peak_day = findall(x->x==maximum(collect(values(daily_p_max))),collect(values(daily_p_max)))[1]   #find day with highest peak production
    highest_avg_p_day = findall(x->x==maximum(collect(values(daily_avg_prod))),collect(values(daily_avg_prod)))[1]   #find day with highest average production

    if daily_avg_prod[highest_avg_p_day] > ref["mean_daily_power_production"]

        ref["mean_daily_power_production"] = daily_avg_prod[highest_avg_p_day]
        ref["peak_power"] = daily_p_max[highest_avg_p_day]
        ref["month"] = month
        ref["day"] = highest_avg_p_day

    end
    #=
    if daily_p_max[highest_peak_day]>ref["peak_power"]  #Alternative condition based on Peak Power
        ref["peak_power"] = daily_p_max[highest_peak_day]
        ref["mean_daily_power_production"] = daily_avg_prod[highest_peak_day]
        ref["month"] = month
        ref["day"] = highest_peak_day
    end 
    =#

end

# Criteria for the selection of the pv generation profile during winter (i.e. least productive)
# Default is based on lowest daily average production, other optional criteria is lowest peak power
function winter_criteria(monthly_profile::OrderedDict, month::Int64, ref::Dict)

    daily_avg_prod, daily_p_max = calc_daily_data(monthly_profile)

    highest_peak_day = findall(x->x==minimum(collect(values(daily_p_max))),collect(values(daily_p_max)))[1]   #find day with highest peak production
    highest_avg_p_day = findall(x->x==minimum(collect(values(daily_avg_prod))),collect(values(daily_avg_prod)))[1]   #find day with highest average production

    if daily_avg_prod[highest_avg_p_day] < ref["mean_daily_power_production"]
        
        ref["mean_daily_power_production"] = daily_avg_prod[highest_avg_p_day]
        ref["peak_power"] = daily_p_max[highest_avg_p_day]
        ref["month"] = month
        ref["day"] = highest_avg_p_day

        
    end
    #=
    if daily_p_max[highest_peak_day] < ref["peak_power"]  #Alternative condition based on Peak Power
        ref["peak_power"] = daily_p_max[highest_peak_day]
        ref["mean_daily_power_production"] = daily_avg_prod[highest_peak_day]
        ref["month"] = month
        ref["day"] = highest_peak_day
    end 
    =#
end

# Compute the monthly average power production of a pv profile 
function calc_avg_monthly_power(daily_profile::OrderedDict)

    p_prod_list = []

    for (day,p_prod) in daily_profile

        [push!(p_prod_list, round(k, digits = 5)) for k in p_prod]
       
    end

    return mean(p_prod_list)
end

# Compute the daily average power productuon and max power peak of a monhtly production profile 
function calc_daily_data(daily_profile::OrderedDict)
    daily_avg_prod = OrderedDict{Int64, Float64}()
    daily_p_max = OrderedDict{Int64, Float64}()

    for (day,p_prod) in daily_profile
        
        daily_avg_prod[day] =  round(mean(p_prod), digits = 5)
        daily_p_max[day] = round(maximum(p_prod), digits = 5)
       
    end

    return daily_avg_prod, daily_p_max

end

# Computes loading of each branch in each multinetwork (adding it to net_data) 
function mn_calc_branch_loading(data::Dict; congestion_limit::Float64=0.0)


    for (nwid,net) in data["nw"]
        [data["nw"][nwid]["branch"][i]["loading"] = round(max( abs(complex(branch["pt"],branch["qt"]))/branch["rate_a"] , abs(complex(branch["pf"],branch["qf"]))/branch["rate_a"]), digits = 3)*100 for (i,branch) in net["branch"]]
    end

    branch_loading = Dict{String,Any}(string(parse(Int64,nwid)-1) => Dict() for nwid in keys(data["nw"]))

    for (nwid,net) in data["nw"]

        hour = string(parse(Int64,nwid)-1)

        for (id,branch) in net["branch"]

            if branch["loading"] >= congestion_limit
            
                push!(branch_loading[hour], id => branch["loading"]) 

            end

        end

    end

    return branch_loading
end

function find_max_branch_loading(branch_loading::Dict)
    
    relevant_dict = OrderedDict{String,Any}()

    for (nwid, branches) in branch_loading

        loadings = collect(values(branches))
        loading_max = maximum(loadings)

        relevant_dict[nwid] = loading_max

    end

    return relevant_dict

end

# Compute the flexibility offered by each node for each stage of the multinetwork
function mn_calc_flexibility_offered_old(data::Dict, result::Dict)

    n = parse.(Int64,keys(data["nw"])).-1

    flex_loads = OrderedDict{String, Dict}(string(h) => Dict() for h in n)
    tot_load = OrderedDict{String, Dict}(string(h) => Dict() for h in n)

    for (nwid, net) in sort(data["nw"])

        hour = string(parse(Int64,nwid)-1)

        p_load = sum([load["load_p"] for (x,load) in result["solution"]["nw"][nwid]["load"]])
        q_load = sum([load["load_q"] for (x,load) in result["solution"]["nw"][nwid]["load"]])

        for (i,load) in result["solution"]["nw"][nwid]["load"]

            p_nominal = net["load"][i]["pd"]
            q_nominal = net["load"][i]["qd"]
    
            bus_id = string(net["load"][i]["load_bus"])
            flex_loads[hour][bus_id] = Dict("diff_real" => 0.0, "flex_p" => 0.0, "diff_imm" => 0.0, "flex_q" => 0.0)
    
            data["nw"][nwid]["bus"][bus_id]["p_flex"] = 0
            data["nw"][nwid]["bus"][bus_id]["q_flex"] = 0
            data["nw"][nwid]["bus"][bus_id]["flex_type"] = "None"
            
            x_p = load["x_p"]
            y_p = load["y_p"]
    
            x_q = load["x_q"]
            y_q = load["y_q"]

            @assert x_p != y_p
    
            if x_p > 10^-6  # upwards active flexibiility 
    
                flex_loads[hour][bus_id]["diff_real"] = x_p
                flex_loads[hour][bus_id]["flex_p"] = round(x_p/p_nominal, digits = 3)*100
    
                data["nw"][nwid]["bus"][bus_id]["p_flex"] = flex_loads[hour][bus_id]["flex_p"]
                data["nw"][nwid]["bus"][bus_id]["flex_type"] = "Upward"
            
            elseif y_p > 10^-6
    
                flex_loads[hour][bus_id]["diff_real"] = - y_p
                flex_loads[hour][bus_id]["flex_p"] = - round(y_p/p_nominal, digits = 3)*100
    
                data["nw"][nwid]["bus"][bus_id]["p_flex"] = flex_loads[hour][bus_id]["flex_p"]
                data["nw"][nwid]["bus"][bus_id]["flex_type"] = "Downward"
            end
    
            if x_q > 10^-6  # upwards reactive flexibiility 
    
                flex_loads[hour][bus_id]["diff_imm"] = x_q
                flex_loads[hour][bus_id]["flex_q"] = round(x_q/q_nominal, digits = 3)*100
    
                data["nw"][nwid]["bus"][bus_id]["q_flex"] = flex_loads[hour][bus_id]["flex_q"]
                data["nw"][nwid]["bus"][bus_id]["flex_type"] = "Upward"
            
            elseif y_q > 10^-6
    
                flex_loads[hour][bus_id]["diff_imm"] = - y_q
                flex_loads[hour][bus_id]["flex_q"] = - round(y_q/q_nominal, digits = 3)*100
    
                data["nw"][nwid]["bus"][bus_id]["q_flex"] = flex_loads[hour][bus_id]["flex_q"]
                data["nw"][nwid]["bus"][bus_id]["flex_type"] = "Downward"
            end
    
            
            v = []
            [push!(v,string(load["load_bus"])) for (ID,load) in data["nw"][nwid]["load"]]
            b = collect(keys(data["nw"][nwid]["bus"]))
            miss = findall(x->x==0, in(v).(b))
            if !isempty(b[miss])
                for idx in b[miss]
    
                    flex_loads[hour][idx] = Dict(
                        "diff_real" => 0.0,
                        "diff_imm" => 0.0,
                        "flex_p" => 0,
                        "flex_q" => 0,
                    )
                    data["nw"][nwid]["bus"][idx]["p_flex"] = 0
                    data["nw"][nwid]["bus"][idx]["q_flex"] = 0
                end
            end
    
        end

        tot_load[hour]["p_demand"] = p_load
        tot_load[hour]["q_demand"] = q_load
    end


    return flex_loads, tot_load

end

# Compute the flexibility offered by each node for each stage of the multinetwork
function mn_calc_flexibility_offered(data::Dict, result::Dict)

    n = parse.(Int64,keys(data["nw"])).-1

    flex_loads = OrderedDict{String, Dict}(string(h) => Dict() for h in n)
    tot_load = OrderedDict{String, Dict}(string(h) => Dict() for h in n)

    for (nwid, net) in sort(data["nw"])

        hour = string(parse(Int64,nwid)-1)

        p_load_actual = sum([load["load_p"] for (x,load) in result["solution"]["nw"][nwid]["load"]])
        q_load_actual = sum([load["load_q"] for (x,load) in result["solution"]["nw"][nwid]["load"]])

        p_load_nominal = sum([load["pd"] for (x,load) in data["nw"][nwid]["load"]])
        q_load_nominal = sum([load["qd"] for (x,load) in data["nw"][nwid]["load"]])

        for (i,load) in result["solution"]["nw"][nwid]["load"]

            p_nominal = net["load"][i]["pd"]
            q_nominal = net["load"][i]["qd"]

            tanφ = tan(acos(p_nominal/sqrt(p_nominal^2+q_nominal^2)))
    
            bus_id = string(net["load"][i]["load_bus"])
            flex_loads[hour][bus_id] = Dict("diff_real" => 0.0, "flex_p" => 0.0, "diff_imm" => 0.0, "flex_q" => 0.0)
    
            data["nw"][nwid]["bus"][bus_id]["p_flex"] = 0
            data["nw"][nwid]["bus"][bus_id]["q_flex"] = 0
            data["nw"][nwid]["bus"][bus_id]["flex_type"] = "None"
            
            x_p = load["x_p"]
            y_p = load["y_p"]
    
            if x_p > 10^-6  # upwards active flexibiility 
    
                flex_loads[hour][bus_id]["diff_real"] = x_p
                flex_loads[hour][bus_id]["flex_p"] = round(x_p/p_nominal, digits = 3)*100

                flex_loads[hour][bus_id]["diff_imm"] = x_p* tanφ
                flex_loads[hour][bus_id]["flex_q"] = round(x_p* tanφ/q_nominal, digits = 3)*100
    
                data["nw"][nwid]["bus"][bus_id]["p_flex"] = flex_loads[hour][bus_id]["flex_p"]
                data["nw"][nwid]["bus"][bus_id]["flex_type"] = "Upward"
                data["nw"][nwid]["bus"][bus_id]["q_flex"] = flex_loads[hour][bus_id]["flex_q"]
                
            
            elseif y_p > 10^-6
    
                flex_loads[hour][bus_id]["diff_real"] = - y_p 
                flex_loads[hour][bus_id]["flex_p"] = - round(y_p/p_nominal, digits = 3)*100

                flex_loads[hour][bus_id]["diff_imm"] = - y_p * tanφ
                flex_loads[hour][bus_id]["flex_q"] = - round(y_p * tanφ/q_nominal, digits = 3)*100
    
                data["nw"][nwid]["bus"][bus_id]["p_flex"] = flex_loads[hour][bus_id]["flex_p"]
                data["nw"][nwid]["bus"][bus_id]["flex_type"] = "Downward"
                data["nw"][nwid]["bus"][bus_id]["q_flex"] = flex_loads[hour][bus_id]["flex_q"]
            end
            
            v = []
            [push!(v,string(load["load_bus"])) for (ID,load) in data["nw"][nwid]["load"]]
            b = collect(keys(data["nw"][nwid]["bus"]))
            miss = findall(x->x==0, in(v).(b))
            if !isempty(b[miss])
                for idx in b[miss]
    
                    flex_loads[hour][idx] = Dict(
                        "diff_real" => 0.0,
                        "diff_imm" => 0.0,
                        "flex_p" => 0,
                        "flex_q" => 0,
                    )
                    data["nw"][nwid]["bus"][idx]["p_flex"] = 0
                    data["nw"][nwid]["bus"][idx]["q_flex"] = 0
                end
            end
    
        end

        tot_load[hour]["p_demand_actual"] = p_load_actual
        tot_load[hour]["q_demand_actual"] = q_load_actual
        tot_load[hour]["p_demand_nominal"] = p_load_nominal
        tot_load[hour]["q_demand_nominal"] = q_load_nominal
    end


    return flex_loads, tot_load

end

# Branch flows in each stage of multinetwork 
function mn_calc_branch_flow(data::Dict)
    pm_data = get_pm_data(data)

    @assert("per_unit" in keys(pm_data) && pm_data["per_unit"])
    @assert(!haskey(pm_data, "conductors"))

    nws = Dict{String,Any}()

    for (i, pm_nw_data) in pm_data["nw"]
        nws[i] = PowerModels._calc_branch_flow_ac(pm_nw_data)
    end

    return Dict{String,Any}(
        "nw" => nws,
        "per_unit" => pm_data["per_unit"],
    )
end

# Power losses in each stage of multinetwork 
function mn_calc_power_losses(data::Dict)

    
    branch_network_losses = Dict{String, Dict{String,Any}}()
    tot_network_losses = Dict{String, Dict{String,Any}}()
 
    for (nwid,net) in data["nw"]
 
        losses = Dict{String,Any}()
        p_loss = []
        q_loss = []

        for (i,branch) in net["branch"]

            if branch["br_status"] != 0
                pt = branch["pt"]
                qt = branch["qt"]
                pf = branch["pf"]
                qf = branch["qf"]

                losses[i] = Dict(
                    "p_loss" => pt+pf,
                    "q_loss" => qt+qf
                )
                push!(p_loss, pt+pf)
                push!(q_loss, qt+qf)
            else
                losses[i] = Dict(
                    "p_loss" => NaN,
                    "q_loss" => NaN
                )
                push!(p_loss, pt+pf)
                push!(q_loss, qt+qf)
            end
 
        end
        push!(branch_network_losses,nwid => Dict("losses" => losses))
        push!(tot_network_losses, string(parse(Int64,nwid)-1) => Dict("p_loss"=>sum(p_loss), "q_loss"=>sum(q_loss)))
         
    end
 
 
    return Dict("nw"=>branch_network_losses, "per_unit" => true), tot_network_losses
 
end

# give feeder_ID with basic data 
function get_feeder_info(data::Dict)
    # INPUT FOR PATH CALC

    down_node, g, map = downstreamcalcs(data)

    extremes = []
    [push!(extremes, i) for (i,j) in down_node if j ==[]]

    tot_paths = [Tuple[]]
    dist = AbstractFloat[]

    len = Dict((map[j["from_b"]],map[j["to_b"]]) => j["length"] for (ind,j) in data["distance"])

    inv_map = inverse_map(data)
    node_dist = []

    tot_paths, dist, node_dist = paths_calc(g, extremes, len, dist, tot_paths, inv_map, file_name, node_dist)

    paths = []
    feeders = []

    mv_busbar = 0

    #COUNT NUMBEER OF FEEDERS
    for (path,line) in enumerate(tot_paths)

        if line[1][2]!=1 && line[1][2] ∉ feeders
            push!(feeders, line[1][2])
        end
    end

    # CREATE DICT CONTAINING ALL DATA RELEVANT FOR FEEDERS
    alphabet = 'A':'Z'
    feeder_ID = Dict( j => Dict("Name" => "Feeder $i", "Paths" => [], "ciao"=>[],"Paths_ID" => [], "Paths_distance"=> [],"Paths_volt" =>[],"vmin" =>[],"vmax" =>[]) for (i,j) in zip('A':alphabet[length(feeders)], feeders) )
    
    # Write tot_paths in a vectorial way. 
    for vec in tot_paths

        feed = []
        for j in 1:length(vec)

            push!(feed,vec[j][1])
            
        end

        push!(paths,feed)  #Rewrite how paths are made

    end

    for (idx,path) in enumerate(paths)

        if path[end]!=1 && path[2] in keys(feeder_ID)
            push!(feeder_ID[path[2]]["Paths"], path)
            push!(feeder_ID[path[2]]["Paths_ID"], idx)
            push!(feeder_ID[path[2]]["Paths_distance"], node_dist[idx])
        end

    end

    # Add ID of branches belonging to the same feeder

    for (ref, f_data) in feeder_ID

        branches = []
    
        for path in f_data["Paths"]
            
            for (ind,j) in enumerate(path)
        
                if ind<length(path)
        
                    t_bus = path[ind]
                    f_bus = path[ind+1]
        
                    for (idx,branch) in data["branch"]
        
                        if f_bus == branch["f_bus"] && t_bus == branch["t_bus"]
                            if idx ∉ branches
                                push!(branches,idx)
                            end
                        end
                    end
                end
            end
        end
    
        feeder_ID[ref]["Branches"] = branches
    
    end

    # Add ID of bus belonging to the same feeder
    for (ref, f_data) in feeder_ID

        nodes = []

        for path in f_data["Paths"]
            [push!(nodes, bus) for bus in path if bus ∉ nodes]
        end

        feeder_ID[ref]["Buses"] = nodes
        mv_busbar = nodes[1]
    end

    return feeder_ID, mv_busbar, paths
end

# Compute & plot voltage profile of each feeder
function calc_voltage_profile(data::Dict, result::Dict, feeder_info::Dict, paths::Vector)
   
    voltage_profile = Dict{String, Dict}()

    for (nwid, net) in data["nw"]

        feeder_ID_1 = copy(feeder_info)
        
        # Volt profile of each path

        path_volt = sort(Dict(i => zeros(length(paths[i])) for i in keys(paths)))

        # GET VECTOR OF VOLTAGE PROFILES

        for (path_ID, buses) in zip(keys(path_volt),paths)

            for (bus, volt) in result["solution"]["nw"][nwid]["bus"]
                
                bus = parse(Int64,bus)
                if bus in buses
                    ind = findall(x-> x==bus, buses)[1]
                    path_volt[path_ID][ind] = volt["vm"]
                end
            end
        end

        for (f_id,f_data) in feeder_ID_1
            for id_of_path in f_data["Paths_ID"]
                push!(feeder_ID_1[f_id]["ciao"], path_volt[id_of_path])
            end
        end


        # Save profiles for each multinetwork
        mock = Dict{Int64,Any}()

        for (id, f_data) in feeder_ID_1
            push!(mock, id => f_data["ciao"])
            feeder_ID_1[id]["ciao"] = []
            feeder_info[id]["ciao"] = []
        end

        push!(voltage_profile, string(parse(Int64,nwid)-1) => mock)

    end

    return voltage_profile
end

# Provide Dict with node providing flexibility for each multinetwork
function mn_calc_flexible_nodes(load::OrderedDict)

    flexible_nodes = OrderedDict{String,Dict{String,Any}}(hour => Dict{String,Any}() for hour in keys(load))

    for (hour,net) in sort(load)

        nwid = string(parse(Int64,hour)+1)
        mock = Dict{String,Dict{String,Float64}}()

        for (bus_id, bus_data) in net
            if bus_data["flex_p"] != 0.0 && bus_data["flex_q"] != 0.0
    
                mock[bus_id] = Dict(
                    "flex_p"=>bus_data["flex_p"], 
                    "flex_q"=>bus_data["flex_q"],
                    "diff_real" =>bus_data["diff_real"],
                    "diff_imm" => bus_data["diff_imm"]
                )
                
                push!(flexible_nodes,hour => mock)
            end
        end
        if get(flexible_nodes, hour, false) == false
            push!(flexible_nodes, hour => Dict())
        end
    end
    return flexible_nodes

end

# Print stuff
function mn_printing_statements(result::Dict, data::Dict, flexible_nodes::OrderedDict, tot_load::OrderedDict, DG_curtailment::Dict)
    println("\n Termination status: ", result["termination_status"])
    print("\n Objective function: ", result["objective"])

    #print status of each multinetwork
    for (nwid, net) in sort(data["nw"])

        println("\n\nNETWORK SOLUTION AT TIME STEP: ", parse(Int64,nwid),"\n\n")
        print_tot_generation(result["solution"]["nw"][nwid], DG_curtailment[nwid])
        print_tot_load(net, tot_load[nwid])
        print_flex_offered(flexible_nodes[nwid])

    end
        
end

# Calc total flexibility offered in each hour
function calc_tot_DR(flexible_nodes::OrderedDict)

    tot_flex = Dict{String, Any}(hour => Dict() for hour in keys(flexible_nodes))

    for (hour, bus) in flexible_nodes
        nwid = string(parse(Int64,hour)+1)
        try
            p_flex_up = round(sum([node["diff_real"] for (nodeid,node) in bus if node["diff_real"]>0]), digits= 5)
            q_flex_up = round(sum([node["diff_imm"] for (nodeid,node) in bus if node["diff_imm"]>0]), digits= 5)

            tot_flex[hour]["p_flex_up"] = p_flex_up
            tot_flex[hour]["q_flex_up"] = q_flex_up

        catch e
            tot_flex[hour]["p_flex_up"] = 0
            tot_flex[hour]["q_flex_up"] = 0
        end
        # Downwards demand flexibility
        try
            p_flex_dwn = round(sum([node["diff_real"] for (nodeid,node) in bus if node["diff_real"]<0]), digits= 5)
            q_flex_dwn = round(sum([node["diff_imm"] for (nodeid,node) in bus if node["diff_imm"]<0]), digits= 5)
            
            tot_flex[hour]["p_flex_dwn"] = p_flex_dwn
            tot_flex[hour]["q_flex_dwn"] = q_flex_dwn
            
        catch e 
            tot_flex[hour]["p_flex_dwn"] = 0
            tot_flex[hour]["q_flex_dwn"] = 0
        end
    end

    return tot_flex
end

function print_flex_offered(flexible_nodes::Dict)

    println("FLEXIBILITY\n")

    if !isempty(flexible_nodes) 
        # Upwards demand flexibility
        try
            p_flex_up = round(sum([node["diff_real"] for (nodeid,node) in flexible_nodes if node["diff_real"]>0]), digits= 5)
            q_flex_up = round(sum([node["diff_imm"] for (nodeid,node) in flexible_nodes if node["diff_imm"]>0]), digits= 5)

            println("Total active upwards power demand flexibility contracted: ", p_flex_up, " MW") 
            println("Total reactive upwards power demand flexibility contracted: ", q_flex_up, " MVar")

        catch e
            println("No upwards flexibility used")
        end


        # Downwards demand flexibility
        try
            p_flex_dwn = round(sum([node["diff_real"] for (nodeid,node) in flexible_nodes if node["diff_real"]<0]), digits= 5)
            q_flex_dwn = round(sum([node["diff_imm"] for (nodeid,node) in flexible_nodes if node["diff_imm"]<0]), digits= 5)
            
            println("Total active downwards power demand flexibility contracted: ", p_flex_dwn, " MW") 
            println("Total reactive downwards power demand flexibility contracted: ", q_flex_dwn, " MVar")
            
        catch e 
            print("No downwards flexibility used")
        end
     
        
    end
end

function print_tot_load(net::Dict, tot_final_load)

    println("LOAD\n")
    p_load_nominal = round(sum([load["pd"] for (load_id,load) in net["load"]]), digits = 5)
    q_load_nominal = round(sum([load["qd"] for (load_id,load) in net["load"]]), digits = 5)

    println("Nominal active demand: ", p_load_nominal," MW")
    println("Nominal active demand: ", q_load_nominal," MVar")
    println()
    println("Final total active demand: ", round(tot_final_load["p_demand"], digits = 5)," MW")
    println("Final total reactive demand: ", round(tot_final_load["q_demand"], digits = 5)," MVar")
    println()
end

function print_tot_generation(result::Dict, dg_curt::Dict)

    println("GENERATION\n")
    println("Active power provided from slack bus: ", round(result["gen"]["1"]["pg"],digits = 5), "MW")
    println("Reactive power provided from slack bus: ", round(result["gen"]["1"]["qg"],digits = 5), "MVar")
    println()
    if length(result["gen"])>1
        DG_p = round(sum([dg["pg"] for (dg_id, dg) in result["gen"] if dg_id != "1"]), digits = 5)
        DG_q = round(sum([dg["qg"] for (dg_id, dg) in result["gen"] if dg_id != "1"]), digits = 5)
        println("Total active power generation from DGs: ",DG_p, " MW")
        println("Total reactive power generation from DGs: ",DG_q, " MVar")

        if !isempty(dg_curt)
            tot_p_curt = round(sum(values["curtail_abs"] for (x,values) in dg_curt), digits = 5)
            println("Total curtailed active power: ", tot_p_curt)
        end
    end

end

# Print hourly statements with multiple representative days 
function hourly_printing_statements(result::Dict, data::Dict, flexible_nodes::OrderedDict, tot_load::OrderedDict, DG_curtailment::Dict, max_branch_loading::OrderedDict, abs_max_min_volt::Dict)

    #print status for each hour
    for (nwid, net) in sort(data["nw"])

        h = string(parse(Int64,nwid)-1) # Nwid in reality is 1 hour in advance (eg nwid = 13 means hour = 12:00)

        println("\n\nNETWORK SOLUTION AT TIME STEP: ", h,"\n\n")
        print_tot_generation(result["solution"]["nw"][nwid], DG_curtailment[h])
        print_tot_load(net, tot_load[h])
        print_flex_offered(flexible_nodes[h])
        print_max_branch_load(max_branch_loading)
        print_max_min_volt(abs_max_min_volt)

    end
        
end

function print_max_branch_load(max_branch_loading::OrderedDict)
    b_max, nw = findmax(collect(values(max_branch_loading)))
    nw = collect(keys(max_branch_loading))[nw]
    println( "Maximum branch loading: $b_max at hour $nw")
end

function print_max_min_volt(abs_max_min_volt::Dict)
    d = abs_max_min_volt
    feeder_v_min = reduce((x, y) -> d[x]["vmin"] ≤ d[y]["vmin"] ? x : y, keys(d))
    feeder_v_max = reduce((x, y) -> d[x]["vmax"] >= d[y]["vmax"] ? x : y, keys(d))
    
    println("Feeder $feeder_v_min: V min = ",d[feeder_v_min]["vmin"]," at ",d[feeder_v_min]["hour_vmin"])
    println("Feeder $feeder_v_max: V max = ",d[feeder_v_min]["vmax"]," at ",d[feeder_v_min]["hour_vmax"])
end

# Struct needed to plot only one label when plotting multiple function

struct onlyone <: AbstractMatrix{Bool}
    v::Bool
end
function Base.iterate(o::onlyone, state=1)
      state == 1 ? o.v : !o.v, state+1
end
Base.size(o::onlyone) = (1,typemax(Int))
Base.length(o::onlyone) = typemax(Int)
Base.getindex(o::onlyone,i) = i == 1 ? o.v : !o.v
Base.getindex(o::onlyone,i,j) = j == 1 ? o.v : !o.v

function pass_vector_data(feeder_ID::Dict, f_id::Int64, n_profile::Vector)

    nrows = maximum(extrema(length, values(feeder_ID[f_id]["Paths_distance"])))   # longest path 
    ncols = length(feeder_ID[f_id]["Paths_distance"])  #number of paths
    y = fill(NaN, nrows, ncols)

    path_volt = copy(n_profile)

    [y[1:length(profile),i] = profile for (i,profile) in enumerate(n_profile)]
        
    return y 
end

function include_v_bounds(condition::Bool)
    if condition
        Plots.plot!([0.9], color=:red ,seriestype = "hline", linewidth = 3, label = "v = 0.9 pu")
        Plots.plot!([1.1], color=:red ,seriestype = "hline", linewidth = 3, label = "v = 1.1 pu")
    end

end

function x_limit_condition(x::Matrix, x_bound::Bool)
    if x_bound
        x0 = x[1,1]  #first value equals 0, corresponding to the HV/MV station
        x1 = x[2,2]  #second value corresponds to the closest node to the primary station
        xmax = maximum(filter(!isnan,x))  # farthest node in the feeder

        ratio = round(x1/xmax, digits = 4)

        println(ratio)
        
        if ratio > 0.4
            xlim_min = 0.4*xmax
            xlim_max = 1.1*xmax
            Plots.plot!(xlims = (xlim_min, xlim_max))  #show plot between 40% and 110% of x range
        end
    end
end


# Plot feeder profile variation within multinetworks 
function mn_plot_united_feeder_voltage_profile(voltage_profile::Dict, feeder_ID::Dict; save::Bool = false, save_path::String="", v_bound::Bool = false, x_bound::Bool = false)
    
    feeder_profile = OrderedDict{Int64,Dict}()

    for (f_id,x) in feeder_ID
        feeder_profile[f_id] = Dict()
        for (nwid, data_stored) in voltage_profile
            push!(feeder_profile[f_id], parse(Int64,nwid) => data_stored[f_id])
        end
    end
    #= 
    hello = get_distance_to_voltage_profile(voltage_profile, feeder_ID)
    for (f_id,x) in feeder_ID
        feeder_profile[f_id] = Dict()
        for (nwid, data_stored) in hello
            push!(feeder_profile[f_id],nwid=> data_stored[f_id])
        end
    end
    =#

    for (f_id, profiles_per_net) in feeder_profile
        
        nrows = maximum(extrema(length, values(feeder_ID[f_id]["Paths_distance"])))   # longest path 
        ncols = length(feeder_ID[f_id]["Paths_distance"])  #number of paths
        x = fill(NaN, nrows, ncols)

        [x[1:length(i),j] = i for (j,i) in enumerate(feeder_ID[f_id]["Paths_distance"])]  #Fill x with the distance values for each path

        y = []
        
        my_colors = get(ColorSchemes.tableau_20, LinRange(0,1,24))  # select colors 
    
        for (nwid, n_profile) in sort(profiles_per_net)  # y is a multimatrix containing the values of the voltage in each multinetwork 
    
            push!(y, pass_vector_data(feeder_ID,f_id, n_profile))
            
        end

        #=
        vmin = minimum(filter(!isnan,y))
        vmax = maximum(filter(!isnan,y))

        feeder_ID[feeder]["vmin"] = vmin
        feeder_ID[feeder]["vmax"] = vmax
        =#
        #label_min = "v_min "* string(round(vmin, digits = 4)) * "pu"
        #label_max = "v_max "* string(round(vmax, digits = 4)) * "pu"

        
        
        plot = Plots.plot(x,y[1]; marker=(:circle,4), color = my_colors[1], linewidth = 2 , label = "h1", primary=onlyone(true))

        for i in 2:length(profiles_per_net) #plot the remaining profiles for multinetworks
        

            Plots.plot!(x,y[i]; marker=(:circle,4), color = my_colors[i], linewidth = 2 , label = "h$i", primary=onlyone(true))
            xlabel!("Distance (km)")
            ylabel!("Voltage drop (p.u.)")
            title!("Voltage drop for: $(feeder_ID[f_id]["Name"])")
            Plots.plot!(size = (1000,800))
            
        end
        
        #=
        plot = Plots.plot(x,y[1]; line_z = y[1], linewidth = 2, label = false)

        for i in 2:length(profiles_per_net) #plot the remaining profiles for multinetworks
        
            Plots.plot!(x,y[i]; line_z = y[i], linewidth = 2, label = false)
            xlabel!("Distance (km)")
            ylabel!("Voltage drop (p.u.)")
            title!("Voltage drop for: $(feeder_ID[f_id]["Name"])")
            Plots.plot!(size = (1000,800))
            
        end
        =#

        x_limit_condition(x, x_bound)
        include_v_bounds(v_bound)

        display(plot)

        if save
            Plots.savefig(save_path*"$(feeder_ID[feeder]["Name"]).png")
        end


    end


end

# Plot the voltage deviation in all feeders along the day. In type you can specify either "boxplot", "dotted area", "errorband"
function mn_plot_voltage_variation(feeder_ID::Dict, voltage_profile::Dict,r_day::Int64, type::String)

    feeder_profile = OrderedDict{Int64,Dict}()
    # Convert the voltage profile as a dict with feeder => hour => voltage profile 
    for (f_id,x) in feeder_ID
        feeder_profile[f_id] = Dict()
        for (nwid, data_stored) in voltage_profile
            push!(feeder_profile[f_id], parse(Int64,nwid) => data_stored[f_id])
        end
    end   

    daily_feeder_profile = Dict{Int64,Dict{Int64,Any}}(f_id => Dict() for f_id in keys(feeder_ID))

    df = DataFrame(voltage_profile =[], feeder_name = [], hour =[])

    # Create a dataframe where for each hour we specify the voltage value of all nodes in each feeder 
    for f_id in keys(feeder_ID)

        for hour in parse.(Int64, keys(voltage_profile))

            a = Dict{String,Float64}()
            # Needed to assign just one voltage magnitude value to the distance 
            for (dist, volt) in zip(reduce(vcat,feeder_ID[f_id]["Paths_distance"]),reduce(vcat,feeder_profile[f_id][hour]))
                a["$dist"] = volt
            end

            push!(daily_feeder_profile[f_id], hour => a)
            append!(df[!, "voltage_profile"], collect(values(a)))
            append!(df[!, "feeder_name"] , reduce(vcat, [feeder_ID[f_id]["Name"] for i in 1:length(collect(values(a)))]))
            append!(df[!, "hour"] , Int.(ones(length(collect(values(a))))*hour))
        end

    end

    if type == "errorband"
        
        df |> 
        @vlplot(x = {
            :hour
        },
        ) +
        @vlplot(
            title = {text = "Voltage variation on day $r_day"},
            mark = {:errorband, extent =:ci},
            width = 1000, height = 600,
            y = {
                :voltage_profile,
                type = "quantitative",
                scale= {zero = false},
                title = "Voltage (p.u.)",
            },
            color = {
                :feeder_name,
                legend = {title = "Feeders profile"},
            },
        ) +
        @vlplot(
            :line,
            title = {text = "Voltage variation on day $r_day"},
            y = {
                "mean(voltage_profile)",
                type = :quantitative,
            },
            color = {
                :feeder_name,
                legend = {title = "Feeders profile"},
            },
        ) |> display

    elseif type == "dotted area"
    
        df |> 
        @vlplot(x = {
            :hour
        },
        ) +
        @vlplot(
            title = {text = "Voltage variation on day $r_day"},
            width = 1000, height = 600,
            mark={:point,opacity=0.3},
            y = {
                "voltage_profile:q",
                scale= {domain = [0.9, 1.1]},
                title = "Voltage (p.u.)",
            },
            color = {
                :feeder_name,
                legend = {title = "Feeders profile"},
            },
        )  + 
        @vlplot(
            :line,
            title = {text = "Voltage variation on day $r_day"},
            y = {
                "mean(voltage_profile)",
                title = "Voltage (p.u.)",
            },
            color = {
                :feeder_name,
                legend = {title = "Feeders profile"},
            },
            
        ) |> display

    elseif type == "boxplot"
        
        df |>
        @vlplot(
            title = {text = "Voltage variation on day $r_day"},
            width = 1000, height = 600,
            mark={:boxplot, extent="min-max"},
            x="hour:o",
            y={
                "voltage_profile:q",
                scale= {domain = [0.9, 1.1]},
                title = "Voltage (p.u.)",
            },
            color = {
                :feeder_name,
                legend = {title = "Feeders profile"},
            },
        ) |> display

    end

end

# These two functions can be used in case you want to get the value of voltage mapped with the distance
function assign_distance_to_voltage(feed_ID::Dict, volt_profile::Dict, nwid::String, f_id::Int64)
   
    ciao = Dict{Float64,Float64}()

    for (path_id, distance) in enumerate(feed_ID[f_id]["Paths_distance"])

        for (idx, value) in enumerate(distance)

            ciao[value] = volt_profile[nwid][f_id][path_id][idx]

        end

    end

    #return Dict{Int64,Dict{Float64,Float64}}(f_id => ciao)
    return sort(ciao)
end

function get_distance_to_voltage_profile(v_profile::Dict, feed_ID::Dict)

    unique_profile = Dict(i => Dict() for i in keys(v_profile))

    for (nwid,x) in v_profile
    
        for f_id in keys(x)
            push!(unique_profile[nwid], f_id => assign_distance_to_voltage(feed_ID,v_profile,nwid, f_id))
        end
    
        
    end

    return unique_profile
    
end

function add_corall_prop(net_data::Dict, nwid::String)

    file_name = net_data["name"]
    file_name  = "Official_"*split(file_name, "_")[2]*".m"

    net_data = net_data["nw"][nwid]
    
    feeder_ID_1, mv_busbar = get_feeder_data(net_data, file_name)
    title  = replace(file_name, "Official_"=>"")
    title = uppercasefirst(replace(title,".m" =>""))

    for (k, feeder) in feeder_ID_1
        for (id, branch) in net_data["branch"]
            if id in feeder["Branches"]
                net_data["branch"][id]["feeder"] = feeder["Name"]
            end
        end
        for (idx, bus) in net_data["bus"]
            if parse(Int64,idx) in feeder["Buses"][2:end]
                net_data["bus"][idx]["feeder"] = feeder["Name"]
            end
        end
    
    end
    
    [net_data["branch"][id]["feeder"] = "HV/MV" for (id,branch) in net_data["branch"] if !haskey(branch,"feeder")]
    [net_data["bus"][id]["feeder"] = "HV/MV" for (id,bus) in net_data["bus"] if !haskey(bus,"feeder")]

    net_data["load"] = Dict()

    return net_data
    
end

#Plot network with branch loading
# node_attribute can be either: "p_flex",  "q_flex", "vm", "flex_type","basic" 
# gen_attribute can be either: "pg", "curtailment", "qg", "basic"
# branch_attribute can be either: "loading", "q_loss", "p_loss", "pt", "basic"

function mn_plot_grid(case::Dict, node_attribute::String, gen_attribute::String, branch_attribute::String, nwid::String; zoom::Bool=false, display_flow::Bool = true, save_fig::Bool = false, save_path::String= "")

    # copy data for modification by plots
    data = deepcopy(case)
    data = add_corall_prop(data, nwid)

    # what components to display
    show_components = ["bus", "branch", "gen"]

    prop_node, prop_gen, prop_br = mn_dict_of_proprieties(data, node_attribute, gen_attribute, branch_attribute)

    plot = powerplot(data, 
                    show_flow = display_flow,
                    flow_arrow_size_range = [650],
                    components = show_components,
                    width = 1000, 
                    height = 1000,

                    branch_data = prop_br[:branch_data],  #branch
                    branch_data_type = prop_br[:branch_data_type],
                    branch_size = prop_br[:branch_size],
                    branch_color = prop_br[:color_range],

                    bus_data = prop_node[:bus_data],  #bus
                    bus_data_type = prop_node[:bus_data_type],
                    bus_size = prop_node[:bus_size],
                    bus_color = prop_node[:color_range],

                    gen_data = prop_gen[:gen_data],  #gen
                    gen_data_type = prop_gen[:gen_data_type],
                    gen_size = prop_gen[:gen_size],
                    gen_color = prop_gen[:color_range]
    )

    # BRANCH
    if !(branch_attribute in ["basic", "corall"])
        #plot.layer[1]["layer"][1]["encoding"]["color"]["field"]="branch_Percent_Loading"
        plot.layer[1]["layer"][1]["encoding"]["color"]["legend"]= Dict("orient"=>"bottom-right")
        plot.layer[1]["layer"][1]["encoding"]["color"]["title"]= prop_br[:title]
        plot.layer[1]["layer"][1]["encoding"]["color"]["scale"]["domain"]= prop_br[:range]
        plot.layer[1]["layer"][1]["encoding"]["color"]["scale"]["range"] = prop_br[:color_range]
        #plot.layer[1]["layer"][1]["encoding"]["size"]=Dict("field"=>"BranchPower", "title"=>"Branch BaseMW", "type"=>"quantitative", "scale"=>Dict("range"=>[3,10]))
    
    end

    # BUS
    if !(node_attribute in ["basic", "flex_type", "corall"])
        plot.layer[3]["encoding"]["color"]["legend"]= Dict("orient"=>"bottom-right", "offset" => -30)
        plot.layer[3]["encoding"]["color"]["title"]= prop_node[:title]
        plot.layer[3]["encoding"]["color"]["scale"]["domain"]= prop_node[:range]
        plot.layer[3]["encoding"]["color"]["scale"]["range"] = prop_node[:color_range]
    end
        
    # GEN
    if !(gen_attribute in ["basic", "corall"])
        #plot.layer[4]["encoding"]["color"]["field"]="gen_Percent_Loading"
        plot.layer[4]["encoding"]["color"]["legend"] = Dict("orient"=>"bottom-right", "offset" => -60)
        plot.layer[4]["encoding"]["color"]["title"] = prop_gen[:title]
        plot.layer[4]["encoding"]["color"]["scale"]["domain"] = prop_gen[:range]
        plot.layer[4]["encoding"]["color"]["scale"]["range"] = prop_gen[:color_range]
        #plot.layer[4]["encoding"]["size"]=Dict("field"=>"GenPower", "title"=>"Gen BaseMW", "type"=>"quantitative", "scale"=>Dict("range"=>[300,1000]))
    end
        
    @set! plot.resolve.scale.size=:independent
    network = uppercasefirst(replace(data["name"],"matpower_"=>""))
    @set! plot.title = network*" Grid Plot"
    
    if branch_attribute == "corall"
        @set! plot.resolve.scale.color=:shared
        plot.layer[1]["layer"][1]["encoding"]["color"]["title"]="Feeders"
        plot.layer[1]["layer"][1]["encoding"]["color"]["legend"]=Dict("clipHeight"=>50, "type" => "symbol", "labelFontSize"=>10, "symbolType" => "circle", "symbolSize" => 1000)
        plot.layer[1]["layer"][1]["encoding"]["color"]["scale"]["range"] = prop_br[:color_range]
    end
    

    if zoom
        PowerPlots.Experimental.add_zoom!(plot)
    end

    if save_fig
        save(save_path*"Grid_plot_mn$nwid.html", plot)
    end

    return plot

end

function corall_checks(data::Dict)
    n_feeders = []
    n_feeders = length([push!(n_feeders,branch["feeder"]) for (id,branch) in data["branch"] if !(branch["feeder"] in n_feeders)])
    
    file_name = data["name"]
    title  = uppercasefirst(replace(file_name, "matpower_"=>""))
    
    if title != "Semiurban"
        color_range = colorscheme2array(ColorSchemes.colorschemes[:tableau_10])[1:n_feeders]
    else
        color_range = colorscheme2array(ColorSchemes.colorschemes[:tableau_20])[1:n_feeders]
    end

    return color_range
end

function mn_dict_of_proprieties(data::Dict, node::String, gen::String, branch::String)
    
    dict_node = Dict{String, Dict}()
    dict_gen = Dict{String, Dict}()
    dict_branch = Dict{String, Dict}()

    bus_size = 90
    branch_size = 3
    gen_size = 170

    f_p_max, f_p_min = dict_find_new(data, "bus", "p_flex")
    f_q_max, f_q_min = dict_find_new(data, "bus", "q_flex")
    v_max, v_min = dict_find_new(data, "bus", "vm")
    curt_max, min = dict_find_new(data, "gen", "curtail_%")
    pg_max, pg_min = dict_find_new(data, "gen", "pg")
    qg_max, qg_min = dict_find_new(data, "gen", "qg")
    b_l_max, b_l_min = dict_find_new(data, "branch", "loading")
    p_loss_max, min = dict_find_new(data, "branch", "p_loss")
    q_loss_max, min= dict_find_new(data, "branch", "q_loss")
    p_f_max = dict_find_new(data, "branch", "pf")

    c_range = corall_checks(data)

    # Attributes for nodes (bus)
    if f_p_min >=0
        dict_node["p_flex"] = Dict{Symbol, Any}(
            :bus_size => bus_size,
            :bus_data =>"p_flex",
            :bus_data_type => "quantitative",
            :range => [f_p_min, f_p_max],
            :title => "Active DR offered (%)",
            :color_range => ["#C0C0C0","#000000"] #grey-black     
        )
    elseif f_p_max <=0
        dict_node["p_flex"] = Dict{Symbol, Any}(
            :bus_size => bus_size,
            :bus_data =>"p_flex",
            :bus_data_type => "quantitative",
            :range => [f_p_min, f_p_max],
            :title => "Active DR offered (%)",
            :color_range => ["#000000","#C0C0C0"] #black - grey        
        )
    elseif f_p_min <0 && f_p_max >0
        dict_node["p_flex"] = Dict{Symbol, Any}(
            :bus_size => bus_size,
            :bus_data =>"p_flex",
            :bus_data_type => "quantitative",
            :range => [f_p_min, f_p_max],
            :title => "Active DR offered (%)",
            :color_range => colorscheme2array(ColorSchemes.colorschemes[:RdYlBu_10]) #red - blue       
        )
    end
    
    if f_q_min >=0
        dict_node["q_flex"] = Dict{Symbol, Any}(
            :bus_size => bus_size,
            :bus_data =>"q_flex",
            :bus_data_type => "quantitative",
            :range => [f_q_min, f_q_max],
            :title => "Reactive DR offered (%)",
            :color_range => ["#C0C0C0","#000000"] #grey-black      
        )
    elseif f_q_max <=0
        dict_node["q_flex"] = Dict{Symbol, Any}(
            :bus_size => bus_size,
            :bus_data =>"q_flex",
            :bus_data_type => "quantitative",
            :range => [f_q_min, f_q_max],
            :title => "Reactive DR offered (%)",
            :color_range => ["#000000","#C0C0C0"] #black - grey        
        )
    elseif f_q_min <0 && f_p_max >0
        dict_node["q_flex"] = Dict{Symbol, Any}(
            :bus_size => bus_size,
            :bus_data =>"q_flex",
            :bus_data_type => "quantitative",
            :range => [f_q_min, f_q_max],
            :title => "Reactive DR offered (%)",
            :color_range => colorscheme2array(ColorSchemes.colorschemes[:PiYG_10]),  #pink - green       
        )
    end
    
    dict_node["vm"] = Dict{Symbol, Any}(
        :bus_size => bus_size,
        :bus_data => "vm",
        :bus_data_type => "quantitative",
        #:range => [v_min,v_max],
        :range => [0.9,1.1],
        :title => "Voltage Magnitude (p.u.)",
        :color_range => ["#FFF5EE","#B0E0E6","#000080"],
    )
    dict_node["flex_type"] = Dict{Symbol, Any}(
        :bus_size => bus_size,
        :bus_data => "flex_type",
        :bus_data_type => "nominal",
    )
    dict_node["basic"] = Dict{Symbol, Any}(
        :bus_size => bus_size,
        :bus_data => "ComponentType",
        :bus_data_type => "nominal",
        :color_range => "#C0C0C0" #"#228b22"  #ForestGreen
    )
    dict_node["corall"] = Dict{Symbol, Any}(
        :bus_size => 15,
        :bus_data => "feeder",
        :bus_data_type => "nominal",
        :color_range => c_range
    )


    # Attributes for gen
    dict_gen["pg"] = Dict{Symbol, Any}(
        :gen_size => gen_size,
        :gen_data => "pg",
        :gen_data_type => "quantitative",
        :range => [pg_min,pg_max],
        :title => "Active power generated (MW)",
        :color_range => ["white","purple"],
    )
    dict_gen["qg"] = Dict{Symbol, Any}(
        :gen_size => gen_size,
        :gen_data => "qg",
        :gen_data_type => "quantitative",
        :range => [qg_min,qg_max],
        :title => "Reactive power generated (MVar)",
        :color_range => ["white","orange"],
    )
    dict_gen["curtailment"] = Dict{Symbol, Any}(
        :gen_size => gen_size,
        :gen_data => "curtail_%",
        :gen_data_type => "quantitative",
        :range => [0,curt_max],
        :title => "DG Curtailment (%)",
        :color_range => colorscheme2array(ColorSchemes.colorschemes[:BuPu_3]),  #white to purple  
    )
    dict_gen["basic"] = Dict{Symbol, Any}(
        :gen_size => gen_size,
        :gen_data => "ComponentType",
        :gen_data_type => "nominal",
        :color_range => "purple"
    )
    dict_gen["corall"] = Dict{Symbol, Any}(
        :gen_size => 100,
        :gen_data => "ComponentType",
        :gen_data_type => "nominal",
        :color_range => "black"
    )

    # Attributes for branch
    dict_branch["loading"] = Dict{Symbol, Any}(
        :branch_size => branch_size,
        :branch_data => "loading",
        :branch_data_type => "quantitative",
        #:range => [0,b_l_max],
        :range => [0,100],
        :title => "Branch loading (%)",
        :color_range => ["green","orange","red"],
    )
    dict_branch["q_loss"] = Dict{Symbol, Any}(
        :branch_size => branch_size,
        :branch_data => "q_loss",
        :branch_data_type => "quantitative",
        :range => [0,q_loss_max],
        :title => "Q loss (MVar)",
        :color_range => ["grey","red"],
    )
    dict_branch["p_loss"] = Dict{Symbol, Any}(
        :branch_size => branch_size,
        :branch_data => "p_loss",
        :branch_data_type => "quantitative",
        :range => [0,p_loss_max],
        :title => "P loss (MW)",
        :color_range => ["grey","red"],
    )
    dict_branch["pf"] = Dict{Symbol, Any}(
        :branch_size => branch_size,
        :branch_data => "pf",
        :branch_data_type => "quantitative",
        :range => [0,p_f_max],
        :title => "Power Flow (MW)",
        :color_range => ["green","red"],
    )
    dict_branch["basic"] = Dict{Symbol, Any}(
        :branch_size => branch_size,
        :branch_data => "ComponentType",
        :branch_data_type => "nominal",
        :color_range => "#87cefa" # lightskyblue #00BFF" #deepSkyBlue
    )
    dict_branch["corall"] = Dict{Symbol, Any}(
        :branch_size => 4,
        :branch_data => "feeder",
        :branch_data_type => "nominal",
        :color_range => c_range
    )

    return dict_node[node],dict_gen[gen], dict_branch[branch]

end

function dict_find_new(d::Dict, category::String, property::String, max::Bool=true, min::Bool=true)

    vector = []
    
    [push!(vector,get(data,property,0)) for (ID, data) in d[category]]

    if max && min
        return maximum(vector), minimum(vector)

    elseif max
        return maximum(vector), nothing

    elseif min
        return nothing, minimum(vector)
        
    end

end

#=
Each generator will get installed a power equal to size_std. 
These generators are added to the model. 
The list of generators are passed as a Dict, where you give a list of buses for each feeders. This is done by calling before get_random_generators
=#

function add_curtailable_generators(net_data::Dict, generators::Dict{Any, Any}, size_std, curt)
    
    x = []
    [append!(x,v) for v in values(generators)]

    for (i,gen) in enumerate(x)
        i+=1
        net_data["gen"]["$i"] = Dict(
            "p_nominal" => size_std, 
            "q_nominal"=>0, 
            "pg" =>size_std, 
            "qg" =>0,
            "pmin" => size_std * (1-curt) ,
            "pmax"=> size_std ,
            "qmin" =>0, 
            "qmax"=>0, 
            "gen_bus" => gen, 
            "gen_status"=>1, 
            "index" => i, 
            "source_id" => ["gen", i]
        )
        net_data["bus"]["$gen"]["bus_type"] = 2
    end
end

# Add generator with curtailment capabilities 
function add_single_curtailable_generator(net_data::Dict, size, gen; curt::Float64=0.0)

    if isa(gen, String)
        gen = parse(Int64, gen)
    end
    
    i = length(net_data["gen"]) + 1

    net_data["gen"]["$i"] = Dict(
        "p_nominal" => size, 
        "q_nominal"=>0, 
        "pg" =>size, 
        "qg" =>0, 
        "pmin" => size * (1-curt), 
        "pmax"=>size, 
        "qmin" =>0, 
        "qmax"=>0, 
        "gen_bus" => gen, 
        "gen_status"=>1, 
        "index" => i, 
        "source_id" => ["gen", i]
    )

    net_data["bus"]["$gen"]["bus_type"] = 2
    
end

# Fin max and min voltages for each feeder
function find_max_min_voltage(volt_prof::Dict)
    
    max_min_volt = Dict{String, Any}()

    for (nwid,net) in volt_prof

        mock = Dict{Int64, Any}(feeder_id => Dict() for feeder_id in keys(net))

        for (f_id, feeder) in net

            v_min = minimum(minimum.(feeder))
            v_max = maximum(maximum.(feeder))

            mock[f_id]["v_min"] = v_min
            mock[f_id]["v_max"] = v_max

        end

        max_min_volt[nwid] = mock

    end

    return max_min_volt
end

# Find curtailment for each gen in each production hour
function mn_calc_curtailment(data::Dict, result::Dict, DG_prod_hours::Vector)

    n = parse.(Int64,keys(data["nw"])).-1 # conversion to int and taking into account hour difference

    feeder_curtailment = Dict(string(hour) => Dict() for hour in n)
    
    @assert !isempty(DG_prod_hours)

    for hour in n[∉(DG_prod_hours).(n)]  #non-generation hours 

        nwid = string(hour+1)  #nwid should be kept as reference 
        hour = string(hour)

        for (gen_id, gen) in result["solution"]["nw"][nwid]["gen"]
        
            push!(feeder_curtailment[hour], gen_id => Dict("curtail_%" => 0 , "curtail_abs" => 0))

        end
    end

    for hour in n[in(DG_prod_hours).(n)]  #Evaluate curtailment only during the hours of production

        nwid = string(hour+1)  #nwid should be kept as reference 
        hour = string(hour)

        for (i, gen) in result["solution"]["nw"][nwid]["gen"]
            if i!="1"
                p_ref = data["nw"][nwid]["gen"][i]["pg_ref"]
                p_final = gen["pg"]
                curtail = p_ref - p_final
        
                if curtail > 10^-4
                    perc_curt = round(curtail/p_ref, digits = 4)*100  # % curtailment

                    push!(feeder_curtailment[hour], i=> Dict("curtail_%" => perc_curt , "curtail_abs" => curtail))
                    data["nw"][nwid]["gen"][i]["curtail_%"] = perc_curt  #update net_data
                    data["nw"][nwid]["gen"][i]["curtail_abs"] = curtail
                    
                else
                    data["nw"][nwid]["gen"][i]["curtail_%"] = 0 
                    data["nw"][nwid]["gen"][i]["curtail_abs"] = 0 
                    push!(feeder_curtailment[hour], i=> Dict("curtail_%" => 0 , "curtail_abs" => 0))
                end
            end
        end
    end

    return feeder_curtailment

end

# Find the maximum and minimu voltage for each feeder happening in the whole simulation period
function find_absolute_max_min_voltage(max_min_volt::Dict)

    reference = Dict{Int64, Dict}(f_id => Dict() for f_id in keys(max_min_volt["0"]))
    vmax = Dict{Int64, Float64}()
    vmin = Dict{Int64, Float64}()

    for f_id in keys(max_min_volt["0"])
        hr_max = reduce((x,y) -> max_min_volt[x][f_id]["v_max"] > max_min_volt[y][f_id]["v_max"] ? x : y, keys(max_min_volt))
        hr_min = reduce((x,y) -> max_min_volt[x][f_id]["v_min"] < max_min_volt[y][f_id]["v_min"] ? x : y, keys(max_min_volt))
        reference[f_id]["vmax"] = round(max_min_volt[hr_max][f_id]["v_max"], digits = 4)
        reference[f_id]["vmin"] = round(max_min_volt[hr_min][f_id]["v_min"], digits = 4)
        reference[f_id]["hour_vmax"] = []
        reference[f_id]["hour_vmin"] = []

    end

    for (hour, volt) in max_min_volt

        for (f_id, data) in volt
            
            if round(data["v_max"], digits = 4) == reference[f_id]["vmax"]
                push!(reference[f_id]["hour_vmax"], hour)
            end  
            if round(data["v_min"], digits = 4) == reference[f_id]["vmin"]
                push!(reference[f_id]["hour_vmin"], hour)
            end
        end
            
    end

    return reference
end

# save hours of production for PV
function find_production_hours(pv_profile::DataFrame)
    return pv_profile[pv_profile[!,2].>.0, :time]
end

function mn_find_production_hours(pv_profiles::Tuple, day::Int64)  # Different version to take into account different production profiles
    prod_h_max = Vector() 

    for profile_df in pv_profiles
        profile = profile_df[!, ["time","PV_Profile_$day"]]
        
        prod_h = find_production_hours(profile)

        if length(prod_h) > length(prod_h_max)
            prod_h_max = find_production_hours(profile)
        end

    end
    return prod_h_max  #return the vector with the longest hours of production 
end

function mn_calc_tot_production(data::Dict, result::Dict)

    DG_prod = Dict{String,Any}()

    for (nwid, sol) in result["solution"]["nw"]

        hour = string(parse(Int64,nwid)-1)

        DG_prod[hour] = state_generation(sol, data["nw"][nwid])

    end
    return DG_prod
end

function state_generation(result::Dict, data::Dict)

    mock = Dict{String,Any}("p_nominal" => Dict(), "p_actual" => Dict())

    if length(result["gen"])>1
        DG_p = round(sum([dg["pg"] for (dg_id, dg) in result["gen"] if dg_id != "1"]), digits = 5)
        DG_q = round(sum([dg["qg"] for (dg_id, dg) in result["gen"] if dg_id != "1"]), digits = 5)

        DG_p_nom = round(sum([dg["pg_ref"] for (dg_id, dg) in data["gen"] if dg_id != "1"]), digits = 5)
        DG_q_nom = round(sum([dg["qg_ref"] for (dg_id, dg) in data["gen"] if dg_id != "1"]), digits = 5)

        mock["p_nominal"] = DG_p_nom
        mock["q_nominal"] = DG_q_nom

        mock["p_actual"] = DG_p
        mock["q_actual"] = DG_q

        return mock

    else
        mock["p_nominal"] = mock["q_nominal"] = mock["p_actual"] = mock["q_actual"] = 0
        return mock
    end
end


function mn_calc_branch_zones(net::Dict, br_loading::Dict; red_lim::Float64 = 80.0, orange_lim::Float64 = 50.0)
    c = Dict( i => Dict("red_zone" => Dict("hours" => [],"tot_congestion_hours" => 0), 
                    "orange_zone" => Dict("hours" => [], "tot_congestion_hours" => 0), 
                    "green_zone" => Dict("hours" => [], "tot_congestion_hours" => 0)
                    )
                    for i in keys(net["branch"])
    )


    for (nwid, data) in br_loading

        for (branch_id, branch_load) in data
        
            if branch_load >= red_lim

                push!(c[branch_id]["red_zone"]["hours"], parse(Int64,nwid))
                c[branch_id]["red_zone"]["tot_congestion_hours"] +=1

            elseif orange_lim <= branch_load < red_lim

                push!(c[branch_id]["orange_zone"]["hours"], parse(Int64,nwid))
                c[branch_id]["orange_zone"]["tot_congestion_hours"] +=1
            else

                push!(c[branch_id]["green_zone"]["hours"], parse(Int64,nwid))
                c[branch_id]["green_zone"]["tot_congestion_hours"] +=1
            end
        end
    end

    return c
end

# PACCO DI FUNZIONI PER PLOTTARE
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

    for prop in prop_to_show #creates df with proprieties to show row by row
        df[idx, prop] = calc_hourly_stuff(PMD, prop)
    end

end

function convert_to_single_df(df::DataFrame, sim_periods::Int)
    power = Vector{Float64}()
    type = []
    for col_name in filter(x -> x!= "hour", names(df))
        power = vcat(power, df[!,col_name])
        type =  vcat(type,reduce(vcat, [col_name for i in 0:sim_periods]))
    end

    hours = repeat(0:sim_periods,length(filter(x -> x!= "hour", names(df))))

    return DataFrame(power = power, type = type, hour = hours)
end

function mn_plot_cumulative_daily_profile(r_day::Int64, sim_periods::Int64, data::Dict, prop_to_show::Vector{String})

    df = DataFrame(hour = 0:sim_periods-1)  #create dataframe where values for each hour will be stored
    proprieties_to_show = Vector{String}()

    for prop in prop_to_show
        prop = aliases(prop)
        push!(proprieties_to_show, prop)
        df[!, prop] = zeros(sim_periods)
    end

    for (nwid, net) in data["nw"]
    
        PMD = PowerModelsDataFrame(net)
        idx = parse(Int64,nwid)

        update_columns(df, PMD, proprieties_to_show, idx)
    end

    push!(df,df[1,:])
    df[end,:hour] = sim_periods
    
    data_df = convert_to_single_df(df, sim_periods)
    
    data_df |> 
    @vlplot(
        :line, 
        title = {text = "Time series day $r_day"},
        color = {
            :type,
            legend = { title = "Line type"}
        },
        width = 1000, height = 600,
        x = {
            :hour,
            title = "Hours"
        },

        y = {
            :power,
            title = "Power (MW)"
        },
    
    ) |> display
    
end



###############################
###### HEAT MAP FUNCTIONS #####
###############################


# Identify max gen min load hour and min gen max load hour 
function critical_hour_max_gen_min_load(data::Dict)
    
    diff_min = 0
    diff_max = 10000
    critical_hour = 0
    
    for (nwid, net) in data["nw"]

        nwid = parse(Int64, nwid)

        if nwid in production_hours

            tot_gen = sum(gen["pg"] for (id,gen) in net["gen"] if id != "1")
            tot_load = sum(load["pd"] for (id,load) in net["load"])

            diff = abs(tot_gen-tot_load)

            if tot_load>tot_gen
                if diff < diff_max
                    
                    diff_max = diff
                    println(diff_max, " ", nwid)
                    critical_hour = nwid

                end
            else
                if diff > diff_min

                    diff_min = diff
                    critical_hour = nwid

                end
            end
            
        end

    end

    return critical_hour

end

function critical_hour_min_gen_max_load(data::Dict)
    
    diff_min = 0
    critical_hour = 0
    
    for (nwid, net) in data["nw"]

        nwid = parse(Int64, nwid)
        tot_load = sum(load["pd"] for (id,load) in net["load"])

        if nwid in production_hours

            tot_gen = sum(gen["pg"] for (id,gen) in net["gen"] if id != "1")

            diff = abs(tot_gen-tot_load)

            if tot_load>tot_gen

                if diff > diff_min
                    
                    diff_min = diff
                    println(diff_min, " ", nwid)
                    critical_hour = nwid

                end
            end
            
        else
            if tot_load > diff_min
                diff_min = tot_load
                critical_hour = nwid
            end
        end

    end

    return critical_hour

end

# Return two dicts, one for zone type: keys = hour, value = number of lines congested
function most_congested_hours(branch_zones::Dict)
    r_hours = Vector{Any}()
    o_hours = Vector{Any}()
    for (branch_id, info) in branch_zones
        if !isempty(info["red_zone"]["hours"])
            push!(r_hours, info["red_zone"]["hours"])
        end
        if !isempty(info["orange_zone"]["hours"])
            push!(o_hours, info["orange_zone"]["hours"])
        end
        
    end
    
    red_hours = DataStructures.counter(reduce(vcat, r_hours))
    orange_hours = DataStructures.counter(reduce(vcat, o_hours))

    return red_hours, orange_hours
end

# value of the most critical hour based on a simple index 
function find_most_congested_hour(r_cong_h, o_cong_h, simulation_periods)
    index = Dict{Int64,Float64}()
    for i in 1:simulation_periods
        v = get(r_cong_h,i,0)*0.7+get(o_cong_h,i,0)*0.7  #index is weighted by the sum of the orange and red congested lines 
        index[i] = v
    end
    return findmax(index)[2]   
end

# Compute the voltage sensitivity of a node/s for a variation in the power injection in that same node 
function calc_volt_ptdf(data::Dict, res_1::Dict, bus::Vector, var::Float64)

    function sensitivity_1(net::Dict, res_1::Dict, bus_id::Int64, var::Float64)  #sensitivity obtained by adding gen

        add_single_generator(net, var, bus_id)

        println(keys(net["gen"]))
        net["per_unit"] = true

        pm = PowerModels.instantiate_model(net, ACPPowerModel, no_flex_with_DG) #build_mn_pf_DR_DGC
        res_2 = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=[])
        update_data!(net, res_2["solution"])

        return (res_2["solution"]["bus"]["$bus_id"]["vm"] - res_1["bus"]["$bus_id"]["vm"])/var
    
    end
    function sensitivity_2(net::Dict, ref::Dict,res_1::Dict, bus_id::Int64, var::Float64)  #sensitivity obtained by modifying load

        bus = ref[:bus_loads][bus_id][1] 
        if var > 0
            net["load"]["$bus"]["pd"]-= var
        else
            net["load"]["$bus"]["pd"]+= var
        end

        net["per_unit"] = true
        res_2 = solve_ac_pf(net, JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
        update_data!(net, res_2["solution"])

        return (res_2["solution"]["bus"]["$bus_id"]["vm"] - res_1["bus"]["$bus_id"]["vm"])/var
    
    end

    ptdf = Dict{Int64,Float64}()
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    for b in bus

        net = deepcopy(data)
        
        #push!(ptdf, b  => sensitivity_1(net, res_1, b,var))
        push!(ptdf, b  => sensitivity_2(net, ref, res_1, b,var))

    end

    return ptdf
end
