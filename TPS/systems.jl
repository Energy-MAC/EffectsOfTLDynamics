cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using TLmodels
using PowerSimulationsDynamics

using PowerSystemsExperiments
const PSE = PowerSystemsExperiments

s = System(joinpath(pwd(), "data/raw_data/WSCC_9bus.raw"))

for l in get_components(PSY.StandardLoad, s)
    transform_load_to_constant_impedance(l)
    l.impedance_active_power = l.impedance_active_power 
    l.impedance_reactive_power = l.impedance_reactive_power 
end

gfm_inj = DynamicInverter(
    "GFM", # stands for "Inverter"
    1.0, # ω_ref,
    PSE.converter_high_power(), #converter
    PSE.VSM_outer_control(), #outer control
    PSE.GFM_inner_control(), #inner control voltage source
    PSE.dc_source_lv(), #dc source
    PSE.pll(), #pll
    PSE.filt(), #filter
)

gfl_inj = DynamicInverter(
    "GFL", # stands for "Inverter"
    1.0, # ω_ref,
    PSE.converter_high_power(), #converter
    PSE.GFL_outer_control(), #outer control
    PSE.GFL_inner_control(), #inner control voltage source
    PSE.dc_source_lv(), #dc source
    PSE.pll(), #pll
    PSE.filt(), #filter
)

sm_inj = DynamicGenerator(
    "SM", # stands for "Generator"
    1.0, # ω_ref,
    PSE.AF_machine(), #machine
    PSE.shaft_no_damping(), #shaft
    PSE.avr_type1(), #avr
    PSE.tg_none(), #tg
    PSE.pss_none(), #pss
)


for g in get_components(Generator, s)
    if get_number(get_bus(g)) == 1
        sm_inj = DynamicGenerator(
        get_name(g), # stands for "Generator"
        1.0, # ω_ref,
        PSE.AF_machine(), #machine
        PSE.shaft_no_damping(), #shaft
        PSE.avr_type1(), #avr
        PSE.tg_none(), #tg
        PSE.pss_none(), #pss
        )
        add_component!(s, sm_inj, g)
    else
        gfm_inj = DynamicInverter(
            get_name(g), # stands for "Inverter"
            1.0, # ω_ref,
            PSE.converter_high_power(), #converter
            PSE.VSM_outer_control(), #outer control
            PSE.GFM_inner_control(), #inner control voltage source
            PSE.dc_source_lv(), #dc source
            PSE.pll(), #pll
            PSE.filt(), #filter
        )
        add_component!(s, gfm_inj, g)
    end
end

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
    get_name(first(get_components(Line, s))),
    10.0,
    1.0,
    1.0
)

sys_statpi = create_statpi_system(s, line_params) 

# to_json(sys_statpi, "sys_statpi", force = true)
sys_dynpi = create_dynpi_system(s, line_params)
sys_MSSB = create_MSSB_system(s, line_params)
sys_MSMB = create_MSMB_system(s, line_params)

for line in get_components(Line, sys_statpi)
    name = line.name
    r = line.r
    x = line.x
    println("Line $name has r = $r, x = $x")
end

for line in get_components(DynamicBranch, sys_dynpi)
    name = get_name(line)
    r = get_r(line)
    x = get_x(line)
    println("Line $name has r = $r, x = $x")
end

for line in get_components(DynamicBranch, sys_MSSB)
    name = get_name(line)
    r = get_r(line)
    x = get_x(line)
    println("Line $name has r = $r, x = $x")
end

for load in get_components(StandardLoad, sys_statpi)
    println(load.bus.name)
    println("P = $(load.impedance_active_power)")
    println("Q = $(load.impedance_reactive_power)")
end

for gen in get_components(Generator, sys_statpi)
    println(gen.bus.name)
    println("P = $(get_active_power(gen))")
    println("Q = $(gen.reactive_power)")
end

using PowerFlows

sol = solve_powerflow(ACPowerFlow(), s)
sol["flow_results"]


for tx in get_components(Transformer2W, sys_statpi)
    println("r = $(tx.r)")
    println("x = $(tx.x)")
end