# Option 1: it solves the connection requests by increasing curtailment progressively. 
# Curt starts at 0. If solution is feasible, then such DG will have allowed_curt = 0. 
# If solution unfeasible, curt is increased by 5%. Loop is stopped for curt = 20% if no feasible solution is found. 

module opt1

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
    using VegaLite, Vega, Query
    using TimerOutputs

    const to = TimerOutput()
    @timeit to "file reading" begin 
        include("C:/Workdir/Develop/PF_simulations/My_functions.jl")
        include("Functions/Mn_functions.jl")
        include("Functions/Conn_functions.jl")

        include("My_ref/My_ref_mn_only_DGC.jl")

        include("Functions/Profiles_functions.jl")
        include("Functions/Solve_conn_requests.jl")

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

        parameters = YAML.load_file("Input.yaml")
        #DG_connection_requests = sort_dict(YAML.load_file("Connection_requests.yaml")["GENERATORS"])
        DG_connection_requests = sort_dict(YAML.load_file("sku.yaml")["GENERATORS"])

    end

    ######################################
    ########### INPUT DATA ###############
    ######################################

    @timeit to "processing input data" begin 
        pv_profiles, pv_profile_key, load_profiles, load_profile_key, frequency_of_occurance, n_representative_days, feeder_ID, paths = 
            elaborate_input_data(net_data,parameters; candidate_gen_nodes = []) 
        #plot_profiles(pv_profiles,load_profiles,pv_profile_key,load_profile_key, frequency_of_occurance)


        AGG_RES = Dict{Int64, Any}()        #Aggregate Results
        AGG_NET_DATA = Dict{Int64, Any}()     #Aggregate net data
        quick_summary = Dict{Int64, OrderedDict}( r_day => OrderedDict() for r_day in keys(frequency_of_occurance))   #Store main takeaways

        requests_status = Dict(DG => true for DG in keys(DG_connection_requests))
    end

    ######################################
    ########## REQUESTS EVALUATION #######
    ######################################


    @timeit to "evaluating connection requests" begin 
        for (DG, DG_data) in DG_connection_requests

            introduce_DG_requests(net_data, DG, DG_data, parameters)  # adding dg request 
            curt = 0
            println("Solving yearly simulation...")
            sol_feasible = false

            while curt <= 20 && sol_feasible == false
                
                sol_feasible = true
                
                update_gen_curt(net_data, DG_data["Connection_node"], curt)

                for repr_day in keys(frequency_of_occurance)

                    # Prepare model, add gen and load profiles
                    mn_net_data, production_hours = prepare_model(
                            net_data, 
                            repr_day, 
                            pv_profiles, 
                            pv_profile_key, 
                            load_profiles, 
                            load_profile_key,
                    )
                
                    # Solve model 
                    mn_net_data, result, flag_inner = solve_model(mn_net_data, build_mn_pf_DGC) 
                    
                    if flag_inner  # if flag_inner=true it means that solution was unfeasible for one repr_day 

                        sol_feasible = false  # solution obtained is unfeasible

                        if curt < 20  # let's try to increase curtailment
                            curt += 5
                            println("Increasing curtailment to $curt %")
                            break
                        else  # curtailment limit reached. DG request cannot be accepted
                            requests_status[DG] = false
                            remove_single_generator(net_data, DG_data)
                            println("\x1b[31mRequest denied. Proceeding to evaluate other requests...\x1b[0m")
                            curt = 1000
                            break
                        end
                    end
                    
                end
            end

            if requests_status[DG]
                DG_connection_requests[DG]["allowed_curt"] = curt/100
                println("\x1b[32mRequest accepted!\x1b[0m")
            end

        end
    end

    println("################################################")
    println("\x1b[32mEVALUATION COMPLETED\x1b[0m")
    println("################################################\n\n")
    println("DGs accepted :")
    [println("- ",DG) for (DG,status) in requests_status if status]

    ######################################
    ####### FINAL CONFIG ANALYISIS #######
    ######################################

    println("Checking final solution:\n\n")

    @timeit to "solving final grid config" begin
        for repr_day in keys(frequency_of_occurance)

            # Prepare model, add gen and load profiles
            mn_net_data, production_hours = prepare_model(
                    net_data, 
                    repr_day, 
                    pv_profiles, 
                    pv_profile_key, 
                    load_profiles, 
                    load_profile_key,
            )

            # Solve model 
            mn_net_data, result, flag_inner = solve_model(mn_net_data, build_mn_pf_DGC) 
            
            # Get informations from results
            flexible_nodes, tot_load, tot_DR, DG_curtailment, DG_production, max_branch_loading, voltage_profile, abs_max_min_volt = evaluate_results(
                                                                                                                    mn_net_data, 
                                                                                                                    result, 
                                                                                                                    production_hours,
                                                                                                                    parameters, 
                                                                                                                    feeder_ID, 
                                                                                                                    paths
            )

            push!(AGG_NET_DATA, repr_day => mn_net_data)
            push!(AGG_RES, repr_day => result)

            #hourly_printing_statements(result, mn_net_data, flexible_nodes, tot_load,DG_curtailment,max_branch_loading, abs_max_min_volt)    
            #mn_plot_cumulative_daily_profile(repr_day, parameters["simulation_periods"], mn_net_data, ["nominal generation", "actual generation", "load nominal"])
            #mn_plot_voltage_variation(feeder_ID, voltage_profile, repr_day, "errorband")
            #mn_plot_united_feeder_voltage_profile(voltage_profile, feeder_ID)
            
            update_summary(quick_summary[repr_day], parameters, result, DG_production, tot_DR, tot_load, abs_max_min_volt, max_branch_loading, frequency_of_occurance[repr_day], production_hours)
        end
    end

    ######################################
    #### IDENTIFYING CRITICAL BRANCHES ###
    ######################################

    @timeit to "branch ranking" begin
        branch_rank, branch_loadings, branch_df = heat_map(net_data, quick_summary, AGG_NET_DATA, parameters)
    end

    #show(to)
end

