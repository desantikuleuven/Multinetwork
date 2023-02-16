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

include("opt1.jl")
include("opt2.jl")

include("C:/Workdir/Develop/PF_simulations/My_functions.jl")
include("Functions/Mn_functions.jl")
include("Functions/Conn_functions.jl")
include("Functions/Profiles_functions.jl")
include("Functions/Solve_conn_requests.jl")

#plot_grid(opt1.net_data,"basic", "basic", "rank";display_flow=false)
#plot_grid(opt2.net_data,"basic", "basic", "rank";display_flow=false)
mn_plot_grid(opt2.AGG_NET_DATA[85],"basic", "curtailment", "loading","12";display_flow = true )

mn_plot_branch_yearly_loading(opt2.quick_summary, opt1.branch_df,"87")