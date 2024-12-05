cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
using InfrastructureSystems
using Plots
using PlotlyJS, DataFrames
using TLmodels
using CSV
using DataFramesMeta
using LaTeXStrings
using Logging
using ZIPE_loads
using Dates

using PowerSystemsExperiments
const PSE = PowerSystemsExperiments

s = System(joinpath(pwd(), "data/raw_data/WSCC_9bus.raw"))

for l in get_components(PSY.StandardLoad, s)
    transform_load_to_constant_impedance(l)
    l.impedance_active_power = l.impedance_active_power 
    l.impedance_reactive_power = l.impedance_reactive_power 
end

gfm_inj() = DynamicInverter(
    "GFM", # stands for "Inverter"
    1.0, # ω_ref,
    PSE.converter_high_power(), #converter
    PSE.VSM_outer_control(), #outer control
    PSE.GFM_inner_control(), #inner control voltage source
    PSE.dc_source_lv(), #dc source
    PSE.pll(), #pll
    PSE.filt(), #filter
)

gfl_inj() = DynamicInverter(
    "GFL", # stands for "Inverter"
    1.0, # ω_ref,
    PSE.converter_high_power(), #converter
    PSE.GFL_outer_control(), #outer control
    PSE.GFL_inner_control(), #inner control voltage source
    PSE.dc_source_lv(), #dc source
    PSE.pll(), #pll
    PSE.filt(), #filter
)

sm_inj() = DynamicGenerator(
    "SM", # stands for "Generator"
    1.0, # ω_ref,
    PSE.AF_machine(), #machine
    PSE.shaft_no_damping(), #shaft
    PSE.avr_type1(), #avr
    PSE.tg_none(), #tg
    PSE.pss_none(), #pss
)

# taken from TLModels.jl `TLmodels_tutorial.ipynb`
impedance_csv = "data/cable_data/dommel_data.csv"
capacitance_csv = "data/cable_data/dommel_data_C.csv"

M = 3
z_km, y_km, z_km_ω, Z_c = get_line_parameters_from_data(impedance_csv, capacitance_csv, M)

line_length_dict = Dict(
    "Bus 5-Bus 4-i_1" => 90,
    "Bus 7-Bus 8-i_1" => 80,
    "Bus 6-Bus 4-i_1" => 100,
    "Bus 7-Bus 5-i_1" => 170,
    "Bus 8-Bus 9-i_1" => 110,
    "Bus 9-Bus 6-i_1" => 180,
)

line_params = LineModelParams(
    z_km, 
    y_km, 
    z_km_ω, 
    Z_c,
    M,
    line_length_dict,    
    "Bus 7-Bus 5-i_1",
    10.0,
    1.0,
    1.0
)
# line_params = LineModelParams(
#     z_km, 
#     y_km, 
#     z_km_ω, 
#     Z_c,
#     M,
#     line_length_dict,    
#     "Bus 6-Bus 4-i_1",
#     10.0,
#     1.0,
#     1.0
# )

function no_change(sys::System, params::LineModelParams)
    return sys
end

function set_power_setpt!(sys::System, scale::Real)
    for load in get_components(StandardLoad, sys)
        set_impedance_active_power!(load, get_impedance_active_power(load)*scale)
        set_current_active_power!(load, get_current_active_power(load)*scale)
        set_constant_active_power!(load, get_constant_active_power(load)*scale)
        
        set_impedance_reactive_power!(load, get_impedance_reactive_power(load)*scale)
        set_current_reactive_power!(load, get_current_reactive_power(load)*scale)
        set_constant_reactive_power!(load, get_constant_reactive_power(load)*scale)
    end
    # if scale <= 1.0 return sys end
    for gen in get_components(Generator, sys)
        if gen.bus.bustype == ACBusTypes.PV
            # set_base_power!(g, g.base_power * p.load_scale)
            set_active_power!(gen, get_active_power(gen) * scale)
            set_reactive_power!(gen, get_reactive_power(gen) * scale)
        end
    end
    return sys
end

### Defining a function to scale the line impedance
# function scale_line_length!(sys::System, scale::Real)
#     for line in get_components(Line, sys)
#         line.r = line.r * scale
#         line.x = line.x * scale
#         line.b = (from = line.b.from*scale, to = line.b.to*scale)
#     end
#     return sys
# end

function small_signal_tripped(gss::GridSearchSys, sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
    if isnothing(sim) return missing end
    sys = deepcopy(sim.sys)
    remove_component!(Line, sys, sim.perturbations[1].branch_name)
    newsim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 1.0), disable_timer_outputs=true)
    sm = small_signal_analysis(newsim)
    return sm
end
                    
cases = [sm_inj() sm_inj() gfl_inj()
        sm_inj() gfl_inj() gfl_inj()
        sm_inj() gfm_inj() gfl_inj()                        
        sm_inj() sm_inj() gfm_inj()
        sm_inj() gfm_inj() gfm_inj()
        gfm_inj() sm_inj() gfm_inj()
        ]

gss = GridSearchSys(s, cases,
                        ["Bus1", "Bus 2", "Bus 3"]) # just make sure the busses are in the right order
set_chunksize!(gss, 200)

line_adders = Dict{String, Function}([
    "statpi"=>create_statpi_system,
    "dynpi"=>create_dynpi_system,
    "MSSB"=>create_MSSB_system,
    "MSMB"=>create_MSMB_system,
])
load_scale_range = collect(0.5:0.5:2.0)
line_scale_range = collect(1.0:0.5:3.0)

line_params_list::Vector{LineModelParams} = []

for scale in line_scale_range
    params_copy = deepcopy(line_params)
    params_copy.line_scale = scale
    push!(line_params_list, params_copy)
end

add_lines_sweep!(gss, line_params_list, line_adders)
add_generic_sweep!(gss, "Load scale", set_power_setpt!, load_scale_range)
# add_generic_sweep!(gss, "Line scale", scale_line_length!, line_scale_range)

# add_result!(gss,
#     ["Bus 1 Injector Current", "Bus 2 Injector Current", "Bus 3 Injector Current"],
#     PSE.get_injector_currents,
# )
# add_result!(gss, "Time", PSE.get_time)
# add_result!(gss, "Eigs", PSE.get_eigenvalues)
# add_result!(gss, "initial_sm", PSE.get_sm)
# add_result!(gss, "final_sm", PSE.small_signal_tripped)
# add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)

current_time_string = string(now())
println(current_time_string)

pathname = "../../../../data/gabrielecr/"*line_params.alg_line_name*"/"*current_time_string
mkdir(pathname)

execute_sims!(gss, BranchTrip(0.5, ACBranch, line_params.alg_line_name), tspan=(0.48, 2.0), dtmax=0.05, run_transient=true, log_path=pathname)