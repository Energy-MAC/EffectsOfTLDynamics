module ExperimentStructs

using PowerSimulationsDynamics
using PowerSystems
using Parameters

@with_kw mutable struct SimParams
    abstol::Float64
    reltol::Float64
    maxiters::Int
    dtmax::Float64
    solver::String
    t_max::Float64
end

@with_kw mutable struct BICParam
    line_name::String
    multiplier::Float64
end

@with_kw mutable struct GenTripParam
    source_type::Type{<:DynamicInjection}
    gen_name::String
end

@with_kw mutable struct CRCParam
    source_type::Type{<:DynamicInjection}
    gen_name::String
    var_to_change::Symbol
    ref_value::Float64
end

@with_kw mutable struct LCParam
    load_type::Type{<:ElectricLoad}
    load_name::String
    var_to_change::Symbol
    ref_value::Float64
end

@with_kw mutable struct LTParam
    load_type::Type{<:ElectricLoad}
    load_name::String
end

@with_kw mutable struct SBVCParam
    var_to_change::Symbol
    ref_value::Float64
end

@with_kw mutable struct PerturbationParams
    t_fault::Float64
    branch_impedance_change_params::Union{BICParam, Nothing} = nothing
    gen_trip_params::Union{GenTripParam, Nothing} = nothing
    crc_params::Union{CRCParam, Nothing} = nothing
    load_change_params::Union{LCParam, Nothing} = nothing
    load_trip_params::Union{LTParam, Nothing} = nothing
    source_bus_voltage_change_params::Union{SBVCParam, Nothing} = nothing
end

@with_kw mutable struct ExpParams
    N::Union{Int64, Nothing}
    M::Union{Int64, Nothing}
    l::Union{Int64, Float64}
    Z_c::Float64
    r_km::Union{Vector{Float64}, Matrix{Float64}}
    x_km::Union{Vector{Float64}, Matrix{Float64}}
    g_km::Union{Vector{Float64}, Matrix{Float64}}
    b_km::Union{Vector{Float64}, Matrix{Float64}}
    r_km_pi::Float64
    x_km_pi::Float64
    Z_c_pi::Float64
    sim_params::SimParams
    perturbation::String
    perturbation_params::PerturbationParams
end

default_bic_params = BICParam("line-102-103", 1.0)
default_gen_trip_params = GenTripParam(DynamicGenerator, "generator-102-1")
default_crc_params = CRCParam(DynamicInverter, "generator-102-1", :V_ref, 0.95)
default_lc_params = LCParam(ElectricLoad, "load-103-1", :P_ref, 1.08)
default_lt_params = LTParam(ElectricLoad, "load-103-1")
default_sbvc_params = SBVCParam(:V_ref, 1.048)

function get_default_perturbation(t_fault::Float64, perturbation::String)
    if perturbation == "BIC"
        perturbation_struct = PerturbationParams(t_fault, default_bic_params, nothing, nothing, nothing, nothing, nothing)
    elseif perturbation == "GenTrip"
        perturbation_struct = PerturbationParams(t_fault, nothing, default_gen_trip_params, nothing, nothing, nothing, nothing)
    elseif perturbation == "CRC"
        perturbation_struct = PerturbationParams(t_fault, nothing, nothing, default_crc_params, nothing, nothing, nothing)
    elseif perturbation == "LoadChange"
        perturbation_struct = PerturbationParams(t_fault, nothing, nothing, nothing, default_lc_params, nothing, nothing)
    elseif perturbation == "LoadTrip"
        perturbation_struct = PerturbationParams(t_fault, nothing, nothing, nothing, nothing, default_lt_params, nothing)
    elseif perturbation == "InfBusChange"
        perturbation_struct = PerturbationParams(t_fault, nothing, nothing, nothing, nothing, nothing, default_sbvc_params)
    else
        return error("Unknown perturbation")
    end
    return perturbation_struct
end

export default_bic_params
export default_gen_trip_params
export default_crc_params
export default_lc_params
export default_lt_params
export default_sbvc_params

export get_default_perturbation

export SimParams
export PerturbationParams
export ExpParams
export BICParam
export GenTripParam
export CRCParam
export LCParam
export LTParam
export SBVCParam

end # module ExperimentStructs
