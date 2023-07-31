using Parameters
using PowerSimulationsDynamics
using PowerSystems

const PSID = PowerSimulationsDynamics
const PSY = PowerSystems

mutable struct SimParams
    abstol::Float64
    reltol::Float64
    maxiters::Int
    dtmax::Float64
    solver::String
    t_max::Float64
end

@with_kw mutable struct PerturbationParams
    t_fault::Float64
    branch_impedance_change_params::BICParam=nothing
    gen_trip::GenTrip=nothing
    crc_params::CRCParam=nothing
    load_change::LCParam=nothing
    load_trip::LTParam=nothing
    source_bus_voltage_change_params::SBVCParam=nothing
end

mutable struct BICParam
    line_name::String
    multiplier::Float64
end

mutable struct GenTrip
    source_type::Type{T<:PSID.DynamicInjection}
    gen_name::String
end

mutable struct CRCParam
    source_type::Type{T<:PSID.DynamicInjection}
    gen_name::String
    var_to_change::Symbol
    ref_value::Float64
end

mutable struct LCParam
    load_type::Type{T<:PSY.ElectricalLoad}
    load_name::String
    var_to_change::Symbol
    ref_value::Float64
end

mutable struct LTParam
    load_type::Type{T<:PSY.ElectricalLoad}
    load_name::String
end

mutable struct SBVCParam
    var_to_change::Symbol
    ref_value::Float64
end

mutable struct ExpParams
    N::Int
    l::Float64
    Z_c::Float64
    r_km::Float64
    x_km::Float64
    g_km::Float64
    b_km::Float64
    sim_params::SimParams
    perturbation::String
    perturbation_params::PerturbationParams
end

export SimParams
export PerturbationParams
export ExpParams
export BICParam
export GenTrip
export CRCParam
export LCParam
export LTParam
export SBVCParam
