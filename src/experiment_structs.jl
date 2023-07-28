using Parameters
using PowerSimulationsDynamics

mutable struct SimParams
    abstol::Float64
    reltol::Float64
    maxiters::Int
    dtmax::FLoat64
    solver::String
    t_max::Float64
end

@with_kw mutable struct PerturbationParams
    t_fault::Float64
    branch_impedance_change_params::BranchImpedanceChange=nothing
    gen_trip::GeneratorTrip=nothing
    crc_params::ControlReferenceChange=nothing
    load_change::LoadChange=nothing
    load_trip::LoadTrip=nothing
    source_bus_voltage_change_params::SourceBusVoltageChange=nothing
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
