using Dates
using CSV
using VegaLite

# Main function to retrieve all profiles
function get_profiles(gen_locations::Vector{String}, load_locations::Vector{String})
    # Gen profiles
    pv_profiles = tuple((get_pv_profiles(pv_location) for pv_location in gen_locations)...)  # Get a tuple of Dicts, each dict is a PV profile in a different location 
    pv_profile_key = Dict(loc => loc_id for (loc_id,loc) in enumerate(gen_locations))

    # Load profiles
    load_profiles = tuple((get_load_profiles(demand_location) for demand_location in load_locations)...)  # Get a tuple of Dicts, each dict is a load profile in a different location 
    load_profile_key = Dict(loc => loc_id for (loc_id,loc) in enumerate(load_locations))

    frequency, days = get_data_info()
    frequency_of_occurance = Dict( d => f for (d,f) in zip(days,frequency))

    return pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance
end

function get_data_info()
    RPF_result_file = "C:/Workdir/Develop/Profiles/Results_2/decision_variables_short.csv" 
    RPF_results = DataFrame(CSV.File(RPF_result_file))
    frequency = RPF_results.weights
    days = RPF_results.periods
    return frequency, days 
end

# Pass all pv profiles to main
function get_pv_profiles(pv_loc::String)

    pv_profile_file_path = "C:/Workdir/Develop/Profiles/PVGIS/PV_profile_2008_"*pv_loc*".json"  # PV generation profile file
    RPF_result_file = "C:/Workdir/Develop/Profiles/Results_2/decision_variables_short.csv"
    RPF_results = DataFrame(CSV.File(RPF_result_file))
    
    return retrieve_pv_profiles(RPF_results, pv_profile_file_path)#, RPF_results.weights
end

function get_daily_pv_profile(ref_day::Int64, pv_production::OrderedDict)
    profile = filter(((k,v),) -> Dates.dayofyear(k) == ref_day,pv_production)  #filter the days of the year based on their number (ie ref_day)
    return collect(values(profile))
end

function retrieve_pv_profiles(RPF_results::DataFrame, file_path::String)
    data = JSON.parsefile(file_path)  # data containing PV profile from PVGIS

    ref_year = data["inputs"]["meteo_data"]["year_max"]
    p_ref = data["inputs"]["pv_module"]["peak_power"]  #kWp

    ref_dates = collect(Dates.DateTime(ref_year):Dates.Hour(1):Dates.DateTime(ref_year+1))
    pop!(ref_dates)
    pv_prod = OrderedDict{DateTime, Float64}(date => (pv["P"]/p_ref)*10^-3 for (date,pv) in zip(ref_dates, data["outputs"]["hourly"]))  #convert data in pu 

    profiles = DataFrame(time = 0:23)

    for day in RPF_results.periods
        profiles[!,"PV_Profile_$day"] = get_daily_pv_profile(day, pv_prod)
    end

    return profiles
end

# Pass all load profiles to main
function get_load_profiles(load_loc::String)

    load_profile_file_path = "C:/Workdir/Develop/Profiles/Load/Load_profile_"*load_loc*".csv"  # PV generation profile file
    RPF_result_file = "C:/Workdir/Develop/Profiles/Results_2/decision_variables_short.csv"
    RPF_results = DataFrame(CSV.File(RPF_result_file))
    
    return retrieve_load_profiles(RPF_results, load_profile_file_path)# RPF_results.weights
end

function get_daily_load_profile(ref_day::Int64, load_profile::DataFrame)
    profile = filter(:Date_Time => x -> Dates.dayofyear(x) == ref_day, load_profile)  #filter the days of the year based on their number (ie ref_day)
    return profile[!,"Global_active_power"]
end

function retrieve_load_profiles(RPF_results::DataFrame, file_path::String)

    load_profile  = CSV.read(file_path, DataFrame, delim = ";")
    load_profile.Date_Time = Dates.DateTime.(load_profile.Date_Time, "dd/mm/yy H:M")
    
    profiles = DataFrame(time = 0:23)
    
    for day in RPF_results.periods
        profiles[!,"Load_Profile_$day"] = get_daily_load_profile(day, load_profile)
    end

    return profiles
end

function assign_gen_pv_profile(data::Dict, locations::Vector{String}; seed::Int64 = 99)
    
    Random.seed!(seed)

    for (gen_id, gen) in data["gen"]
        if gen_id!= "1"
            data["gen"][gen_id]["gen_profile"] =  "PV_profile_"*sample(locations,1)[1] 
        end
    end

end

function assign_load_profile(data::Dict, locations::Vector{String}; seed::Int64 = 99)

    Random.seed!(seed)

    for (load_id, load) in data["load"]
        
        data["load"][load_id]["load_profile"] =  "Load_profile_"*sample(locations,1)[1] 
        
    end

end

# Assign the load profile to each load in the multinetwork network 
function attribute_load_profile_values(data::Dict, load_prof::Tuple, load_profile_key::Dict, r_day::Int64)

    for (n,net) in data["nw"]
        for (i,load) in net["load"]
            h = parse(Int64,n)
            key = replace(load["load_profile"], "Load_profile_" => "")
            idx = load_profile_key[key]
            data["nw"][n]["load"][i]["pd"] = load["pd"] * load_prof[idx][!, "Load_Profile_$r_day"][h]
            data["nw"][n]["load"][i]["qd"] = load["qd"] * load_prof[idx][!, "Load_Profile_$r_day"][h]

            data["nw"][n]["load"][i]["pd_ref"] = data["nw"][n]["load"][i]["pd"]  # used as reference value for curtailment evaluation
            data["nw"][n]["load"][i]["qd_ref"] = data["nw"][n]["load"][i]["qd"]
        end
    end


end

# Assign the load profile to each load in the multinetwork network 
function attribute_gen_profile_values(data::Dict, gen_prof::Tuple, pv_profile_key::Dict, r_day::Int64)

    for (n,net) in data["nw"]
        for (i,gen) in net["gen"]
            if i!= "1"
                h = parse(Int64,n)
                key = replace(gen["gen_profile"], "PV_profile_" => "")
                idx = pv_profile_key[key]
                data["nw"][n]["gen"][i]["pg"] = gen["p_nominal"] * gen_prof[idx][!,"PV_Profile_$r_day"][h]
                data["nw"][n]["gen"][i]["qg"] = gen["p_nominal"] * gen_prof[idx][!,"PV_Profile_$r_day"][h]

                data["nw"][n]["gen"][i]["pg_ref"] = data["nw"][n]["gen"][i]["pg"]  # used as reference value for curtailment evaluation
                data["nw"][n]["gen"][i]["qg_ref"] = data["nw"][n]["gen"][i]["qg"]
            end
        end
    end

end



# PLOT PROFILES

function plot_profiles(pv_profiles::Tuple, load_profiles::Tuple, pv_profile_key::Dict, load_profile_key::Dict, frequency_of_occurance::Dict, aggregate::Bool=false, save_plot::Bool = false; show_load::Bool = false, show_gen::Bool = false, show_all::Bool = false)

    df = collect_data_into_dataframe(pv_profiles, load_profiles, frequency_of_occurance, pv_profile_key, load_profile_key)

    if show_load
        plot_load_profile(df, aggregate, save_plot)
    end  
    if show_gen
        plot_gen_profile(df, aggregate, save_plot)
    end    
    if show_all
        plot_gen_profile(df, aggregate, save_plot)+plot_load_profile(df, aggregate, save_plot)
    end


end

# Gather all data into one big dataframe
function collect_data_into_dataframe(pv_profiles::Tuple, load_profiles::Tuple, fr_of_occurance::Dict, pv_profile_key::Dict, load_profile_key::Dict)

    power = Vector{Float64}()
    type = Vector{String}()

    #collect production data
    for repr_day in keys(fr_of_occurance)
        # Add column of gen profile and its type
        p = reduce(vcat,[pv_profiles[pv_profile_key[i]][!,"PV_Profile_"*string(repr_day)] for i in keys(pv_profile_key)])
        append!(power,p)
        t = reduce(vcat,[reduce(vcat,[loca*"_PV_day_"*string(repr_day) for i in 0:23]) for loca in keys(pv_profile_key)])
        append!(type,t)
        # Add column of gen profile and its type (to be changed when we'll have multiple load profiles WITH SAME SYNTAX)
        p = reduce(vcat,[load_profiles[load_profile_key[i]][!,"Load_Profile_"*string(repr_day)] for i in keys(load_profile_key)])
        append!(power,p)
        t = reduce(vcat,[reduce(vcat,[loca*"_Load_day_"*string(repr_day) for i in 0:23]) for loca in keys(load_profile_key)])
        append!(type,t)

    end

    repetitions = (length(keys(pv_profile_key))+length(keys(load_profile_key)))*n_representative_days  # (NUMBER GEN OF PROFILES + NUMBER OF LOAD PROFILES)*N_REPR_DAYS
    time = reduce(vcat,[[i for i in 0:23] for x in 1:repetitions])
    df = DataFrame(time = time, power = power, type = type)

    return df
end

# PLOT LOAD PROFILES

function plot_load_profile_aggregate(data::DataFrame, save_plot::Bool)

    data.location = popfirst!.(split.(data.type,"_"))

    data |> @vlplot(
        title={text="Average load profiles"},
        :line, 
        color = {
            :location,
            legend = { title = "Average load profile"}
        },
        width = 1000, height = 600,
        x = {
            :time,
            title = "Hours"
        },

        y = {
            "average(power)",
            title = "Power (pu)"
        },
            
    ) |> display
end

function plot_load_profile(data::DataFrame, aggregate::Bool, save_plot::Bool)
    load_data = filter(:type => x -> split(x,"_")[2] == "Load", data)  

    if aggregate
        plot_load_profile_aggregate(load_data, save_plot)
    else

        p = load_data |> @vlplot(
            title={text="Load profiles"},
            :line, 
            color = {
                :type,
                legend = { title = "Line type"}
            },
            width = 1000, height = 600,
            x = {
                :time,
                title = "Hours"
            },
    
            y = {
                "power:q",
                title = "Power (pu)"
            },
            
        ) |> display 
        

    end
end

# PLOT GEN PROFILES

function plot_gen_profile_aggregate(data::DataFrame, save_plot::Bool)
    
    data.location = popfirst!.(split.(data.type,"_"))

    data |> @vlplot(
        title={text="Average generation profiles"},
        :line, 
        color = {
            :location,
            legend = { title = "Average load profile"}
        },
        width = 1000, height = 600,
        x = {
            :time,
            title = "Hours"
        },

        y = {
            "average(power)",
            title = "Power (pu)"
        },
            
    ) |> display
end

function plot_gen_profile(data::DataFrame, aggregate::Bool, save_plot::Bool)
    
    gen_data = filter(:type => x -> split(x,"_")[2] == "PV", data) # In the second argument you can either have PV or Load 

    if aggregate
        plot_gen_profile_aggregate(gen_data, save_plot)
    else

        gen_data |> @vlplot(
            title={text="Generation profiles"},
            :line, 
            color = {
                :type,
                legend = { title = "Line type"}
            },
            width = 1000, height = 600,
            x = {
                :time,
                title = "Hours"
            },

            y = {
                "power:q",
                title = "Power (pu)"
            },
            
        ) |> display

    end
end


# Print stuff
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
    


# Calc total flexibility offered in each hour
function calc_tot_DR(flexible_nodes::OrderedDict)

    tot_flex = Dict{String, Any}(i=> Dict() for i in keys(flexible_nodes))

    for (nwid, bus) in flexible_nodes
        try
            p_flex_up = round(sum([node["diff_real"] for (nodeid,node) in bus if node["diff_real"]>0]), digits= 5)
            q_flex_up = round(sum([node["diff_imm"] for (nodeid,node) in bus if node["diff_imm"]>0]), digits= 5)

            tot_flex[nwid]["p_flex_up"] = p_flex_up
            tot_flex[nwid]["q_flex_up"] = q_flex_up

        catch e
            tot_flex[nwid]["p_flex_up"] = 0
            tot_flex[nwid]["q_flex_up"] = 0
        end
        # Downwards demand flexibility
        try
            p_flex_dwn = round(sum([node["diff_real"] for (nodeid,node) in bus if node["diff_real"]<0]), digits= 5)
            q_flex_dwn = round(sum([node["diff_imm"] for (nodeid,node) in bus if node["diff_imm"]<0]), digits= 5)
            
            tot_flex[nwid]["p_flex_dwn"] = p_flex_dwn
            tot_flex[nwid]["q_flex_dwn"] = q_flex_dwn
            
        catch e 
            tot_flex[nwid]["p_flex_dwn"] = 0
            tot_flex[nwid]["q_flex_dwn"] = 0
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

function update_summary(summary::Dict, result::Dict, mn_net_data::Dict, tot_DR::Dict)

    summary["OF"] = result["objective"]
    summary["Curtailment"] = sum(sum(get(gen,"curtail_abs",0) for (gen_id,gen) in data["gen"]) for (nwid,data) in mn_net_data["nw"])  #total curtailment
    summary["Up_flex"] = sum(flex["p_flex_up"] for (hour,flex) in tot_DR)
    summary["Down_flex"] = sum(flex["p_flex_dwn"] for (hour,flex) in tot_DR)

end

