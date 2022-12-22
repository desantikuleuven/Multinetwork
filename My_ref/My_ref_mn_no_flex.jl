using Memento
const _LOGGER = Memento.getlogger(@__MODULE__)
var(pm::AbstractPowerModel, nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, nw)
var(pm::AbstractPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pm_it_sym, nw, key)
var(pm::AbstractPowerModel, nw::Int, key::Symbol, idx) = _IM.var(pm, pm_it_sym, nw, key, idx)
var(pm::AbstractPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key; nw = nw)
var(pm::AbstractPowerModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key, idx; nw = nw)

#=

    This method is used to compute AC PF with NO flexibility, having DGs as PQ buses

=#


function no_flex_with_DG(pm::AbstractPowerModel)
    for (n,network) in nws(pm)

        variable_bus_voltage(pm, nw=n, bounded = false)    # VOLTAGE CONSTRAINT ACTIVE
        variable_gen_power(pm, nw=n, bounded = false)
        variable_branch_power(pm, nw=n, bounded = false)

        constraint_model_voltage(pm, nw=n)

        for (i,bus) in ref(pm, :ref_buses,nw=n)
            @assert bus["bus_type"] == 3

            constraint_theta_ref(pm, i,nw=n)
            constraint_voltage_magnitude_setpoint(pm, i,nw=n)

        end

        for (i,bus) in ref(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)

            # PQ Bus Constraints
            if length(ref(pm, :bus_gens, i,nw=n)) > 0 && !(i in ids(pm,:ref_buses,nw=n))
                # this assumes inactive generators are filtered out of bus_gens
                @assert bus["bus_type"] == 2  # Should be ==1 but PowerModels is stupid

                for j in ref(pm, :bus_gens, i, nw=n)
                    constraint_gen_setpoint_active(pm, j ,nw=n)
                    constraint_gen_setpoint_reactive(pm, j ,nw=n)
                end
            end
        end

        for i in ids(pm, :branch,nw=n)
            constraint_ohms_yt_to(pm,i,nw=n)
            constraint_ohms_yt_from(pm,i,nw=n)

            constraint_thermal_limit_from_me(pm, i,nw=n)
            constraint_thermal_limit_to_me(pm, i,nw=n)
        end
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