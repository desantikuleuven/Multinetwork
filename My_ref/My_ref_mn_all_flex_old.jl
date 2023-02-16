
using Memento
const _LOGGER = Memento.getlogger(@__MODULE__)
var(pm::AbstractPowerModel, nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, nw)
var(pm::AbstractPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pm_it_sym, nw, key)
var(pm::AbstractPowerModel, nw::Int, key::Symbol, idx) = _IM.var(pm, pm_it_sym, nw, key, idx)
var(pm::AbstractPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key; nw = nw)
var(pm::AbstractPowerModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key, idx; nw = nw)



#=
INFO:

Multinetork model exploiting DG curtailment and DR. ALL CONSTRAINTS INCLUDED.

=#

function build_mn_pf_DR_DGC(pm::AbstractPowerModel)

    for (n,network) in nws(pm)

        variable_bus_voltage(pm, nw=n, bounded = false)
        variable_gen_power_me(pm, nw=n, bounded = true)
        variable_branch_power(pm, nw=n, bounded = false)
        variable_load(pm, nw=n, bounded=true)
        variable_dummy(pm, nw=n)

        constraint_model_voltage(pm, nw=n)  #do nothing

        for (i,bus) in ref(pm, :ref_buses, nw=n)
            @assert bus["bus_type"] == 3

            constraint_theta_ref(pm, i, nw=n)
            constraint_voltage_magnitude_setpoint(pm, i, nw=n)

        end

        for (i,bus) in ref(pm, :bus, nw=n)
            constraint_power_balance_me(pm, i, nw=n)

            constraint_voltage_magnitude_lower_me(pm, i, nw=n)
            constraint_voltage_magnitude_upper_me(pm, i, nw=n)

        end

        for i in ids(pm, :branch, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_thermal_limit_from_me(pm, i, nw=n)
            constraint_thermal_limit_to_me(pm, i, nw=n)
        end

        for i in ids(pm, :load, nw=n)
            constraint_load_factor(pm, i, nw=n)
            constraint_dummy(pm, i, nw=n) 
        end

    end

    objective_min_load_gen_variations(pm)

end



function variable_gen_power_me(pm::AbstractPowerModel; kwargs...)
    variable_gen_power_real_me(pm; kwargs...)
    variable_gen_power_imaginary_me(pm; kwargs...)
end

function variable_gen_power_real_me(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pg = var(pm, nw)[:pg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg")
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if i!=1
                JuMP.set_lower_bound(pg[i], (1-gen["allowed_curt"]) * gen["pg_ref"])
                JuMP.set_upper_bound(pg[i], gen["pg_ref"])
            end
        end
    end

    report && sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end

function variable_gen_power_imaginary_me(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qg = var(pm, nw)[:qg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qg",
        start = comp_start_value(ref(pm, nw, :gen, i), "qg")
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if i!=1
                JuMP.set_lower_bound(qg[i], (1-gen["allowed_curt"]) * gen["qg_ref"])
                JuMP.set_upper_bound(qg[i], gen["qg_ref"])
            end
        end
    end

    report && sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end

function variable_load(pm::AbstractACPModel; kwargs...)
    variable_load_real(pm; kwargs...)
    variable_load_imm(pm; kwargs...)
end

function variable_load_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    load_p= var(pm,nw)[:load_p] = JuMP.@variable(pm.model, 
    [i in ids(pm, nw, :load)],
    base_name="$(nw)_load_p",  start = comp_start_value(ref(pm, nw, :load, i), "pd"))

    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(load_p[i], (1-load["flex_%"]) * load["pd_ref"]) # felxibility 
            JuMP.set_upper_bound(load_p[i], (1+load["flex_%"]) * load["pd_ref"])
        end
    end

    report && sol_component_value(pm, nw, :load, :load_p, ids(pm, nw, :load), load_p)

end

function variable_load_imm(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    
    load_q = var(pm,nw)[:load_q] = JuMP.@variable(pm.model, [i in ids(pm, nw, :load)],
    base_name="$(nw)_load_q",  start = comp_start_value(ref(pm, nw, :load, i), "qd"))

    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(load_q[i], (1-load["flex_%"]) * load["qd_ref"]) # felxibility 
            JuMP.set_upper_bound(load_q[i], (1+load["flex_%"]) * load["qd_ref"])
        end
    end

    report && sol_component_value(pm, nw, :load, :load_q, ids(pm, nw, :load), load_q)

end

function variable_dummy(pm::AbstractPowerModel; nw::Int=nw_id_default, report::Bool=true)

    y_p = var(pm,nw)[:y_p] = JuMP.@variable(pm.model, [i in ids(pm, nw, :load)])
    x_p = var(pm,nw)[:x_p] = JuMP.@variable(pm.model, [i in ids(pm, nw, :load)])
    y_q = var(pm,nw)[:y_q] = JuMP.@variable(pm.model, [i in ids(pm, nw, :load)])
    x_q = var(pm,nw)[:x_q] = JuMP.@variable(pm.model, [i in ids(pm, nw, :load)])

    for (i, load) in ref(pm, nw, :load)
        JuMP.set_lower_bound(y_p[i], 0) # felxibility 
        JuMP.set_lower_bound(x_p[i], 0)
        JuMP.set_lower_bound(y_q[i], 0)
        JuMP.set_lower_bound(x_q[i], 0)
    end

    report && sol_component_value(pm, nw, :load, :x_p, ids(pm, nw, :load), x_p)
    report && sol_component_value(pm, nw, :load, :x_q, ids(pm, nw, :load), x_q)
    report && sol_component_value(pm, nw, :load, :y_p, ids(pm, nw, :load), y_p)
    report && sol_component_value(pm, nw, :load, :y_q, ids(pm, nw, :load), y_q)

end




function objective_min_load_gen_variations(pm::AbstractPowerModel; nw = nw_id_default, kwargs...)
    
    return JuMP.@objective(pm.model, Min, 
        sum(
            sum(gen["pg"]-var(pm,nw, :pg, i) for (i,gen) in ref(pm, nw, :gen) if i!=1) +
            sum(gen["qg"]-var(pm,nw, :qg, i) for (i,gen) in ref(pm, nw, :gen) if i!=1) +
            sum(var(pm,nw, :x_p, i) + var(pm,nw, :y_p, i) for i in ids(pm, nw, :load)) + 
            sum(var(pm,nw, :x_q, i) + var(pm,nw, :y_q, i) for i in ids(pm, nw, :load))
            for (nw, network) in nws(pm))
    )
end



function constraint_power_balance_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
   
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance_me_2(pm, nw, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
end

function constraint_power_balance_me_2(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
    vm   = var(pm, n, :vm, i)
    p    = get(var(pm, n),    :p, Dict()); PowerModels._check_var_keys(p, bus_arcs, "active power", "branch") #real power flowing from all branches (i to j)
    q    = get(var(pm, n),    :q, Dict()); PowerModels._check_var_keys(q, bus_arcs, "reactive power", "branch")
    
    load_p = get(var(pm, n), :load_p, Dict())
    load_q = get(var(pm, n), :load_q, Dict())
    
    pg   = get(var(pm, n),   :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")  #real power injected by generator
    qg   = get(var(pm, n),   :qg, Dict()); PowerModels._check_var_keys(qg, bus_gens, "reactive power", "generator")


    # the check "typeof(p[arc]) <: JuMP.NonlinearExpression" is required for the
    # case when p/q are nonlinear expressions instead of decision variables
    # once NLExpressions are first order in JuMP it should be possible to
    # remove this.
    nl_form = length(bus_arcs) > 0 && (typeof(p[iterate(bus_arcs)[1]]) <: JuMP.NonlinearExpression)

    if !nl_form
        cstr_p = JuMP.@constraint(pm.model,
            sum(p[a] for a in bus_arcs)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(load_p[l] for l in bus_loads)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    else
        cstr_p = JuMP.@NLconstraint(pm.model,
            sum(p[a] for a in bus_arcs)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(load_p[l] for l in bus_loads)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    end

    if !nl_form
        cstr_q = JuMP.@constraint(pm.model,
            sum(q[a] for a in bus_arcs)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(load_q[l] for l in bus_loads)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    else
        cstr_q = JuMP.@NLconstraint(pm.model,
            sum(q[a] for a in bus_arcs)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(load_q[l] for l in bus_loads)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    end

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

function constraint_thermal_limit_from_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        p_fr = var(pm, nw, :p, f_idx)
        q_fr = var(pm, nw, :q, f_idx)
        
        JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (branch["cong_cap"]*branch["rate_a"])^2)  #80% loading 
    end
end

function constraint_thermal_limit_to_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        p_to = var(pm, nw, :p, t_idx)
        q_to = var(pm, nw, :q, t_idx)

        JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (branch["cong_cap"] * branch["rate_a"])^2)  #80% loading 
    end
end

function constraint_voltage_magnitude_upper_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    v_max = bus["vmax"]
    #v_max = 1.2

    v = var(pm, nw, :vm, i)
    JuMP.@constraint(pm.model, v <= v_max)

end

function constraint_voltage_magnitude_lower_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    v_min = bus["vmin"]
    #v_min = 0.8
   

    v = var(pm, nw, :vm, i)
    JuMP.@constraint(pm.model, v_min <= v)

end

function constraint_load_factor(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    load = ref(pm,nw,:load,i)
    load_p = get(var(pm, nw), :load_p, Dict())
    load_q = get(var(pm, nw), :load_q, Dict())

    cosφ = load["pd"]/sqrt(load["pd"]^2+load["qd"]^2)
    φ = acos(cosφ)
    tanφ = tan(φ)

    pd = load_p[i]
    qd = load_q[i]

    @constraint(pm.model, pd*tanφ == qd)
end

function constraint_dummy(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    load = ref(pm,nw,:load,i)
    load_p = get(var(pm, nw), :load_p, Dict())
    load_q = get(var(pm, nw), :load_q, Dict())
    y_p = get(var(pm,nw), :y_p, Dict())
    x_p = get(var(pm,nw), :x_p, Dict())
    y_q = get(var(pm,nw), :y_q, Dict())
    x_q = get(var(pm,nw), :x_q, Dict())

    p = load["pd"]
    p_flex = load_p[i]
    JuMP.@constraint(pm.model, p - p_flex <= y_p[i])
    JuMP.@constraint(pm.model, p_flex - p <= x_p[i])

    q = load["qd"]
    q_flex = load_q[i]
    JuMP.@constraint(pm.model, q - q_flex <= y_q[i])
    JuMP.@constraint(pm.model, q_flex - q <= x_q[i])
end

