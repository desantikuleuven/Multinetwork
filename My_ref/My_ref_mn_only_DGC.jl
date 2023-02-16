
using Memento
const _LOGGER = Memento.getlogger(@__MODULE__)
var(pm::AbstractPowerModel, nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, nw)
var(pm::AbstractPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pm_it_sym, nw, key)
var(pm::AbstractPowerModel, nw::Int, key::Symbol, idx) = _IM.var(pm, pm_it_sym, nw, key, idx)
var(pm::AbstractPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key; nw = nw)
var(pm::AbstractPowerModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key, idx; nw = nw)



#=
INFO:

Multinetork model exploiting only DG curtailment. ALL CONSTRAINTS INCLUDED.

=#

function build_mn_pf_DGC(pm::AbstractPowerModel)

    for (n,network) in nws(pm)

        variable_bus_voltage(pm, nw=n, bounded = false)
        variable_gen_power_me(pm, nw=n, bounded = true)
        variable_branch_power(pm, nw=n, bounded = false)

        constraint_model_voltage(pm, nw=n)  #do nothing

        for (i,bus) in ref(pm, :ref_buses, nw=n)
            @assert bus["bus_type"] == 3

            constraint_theta_ref(pm, i, nw=n)
            constraint_voltage_magnitude_setpoint(pm, i, nw=n)

        end

        for (i,bus) in ref(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)

            constraint_voltage_magnitude_lower_me(pm, i, nw=n)
            constraint_voltage_magnitude_upper_me(pm, i, nw=n)

        end

        for i in ids(pm, :branch, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_thermal_limit_from_me(pm, i, nw=n)
            constraint_thermal_limit_to_me(pm, i, nw=n)
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


function objective_min_load_gen_variations(pm::AbstractPowerModel; nw = nw_id_default, kwargs...)
    
    return JuMP.@objective(pm.model, Min, 
        sum(
            sum(gen["pg"]-var(pm,nw, :pg, i) for (i,gen) in ref(pm, nw, :gen) if i!=1) for (nw, network) in nws(pm)
        )
    )
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

    v = var(pm, nw, :vm, i)
    JuMP.@constraint(pm.model, v <= v_max)

end

function constraint_voltage_magnitude_lower_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    v_min = bus["vmin"]

    v = var(pm, nw, :vm, i)
    JuMP.@constraint(pm.model, v_min <= v)

end

