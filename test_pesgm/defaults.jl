using PowerSimulationsDynamics
using PowerSystems
using Parameters
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics;

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 0.1,
)

line_dict = default_2_bus_line_dict

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

# Perturbation type not relevant for stability analysis 
perturbation_type = "CRC";
t_fault = 0.25;
perturbation_params = get_default_perturbation(t_fault, perturbation_type);

### Extract line data from files for M=1 and M=3
M = 1;
factor_z = 1.0; 
factor_y = 1.0;
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y);

M = 3;
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y);

# Calculate SIL 
V_nom = 230 # kV
Z_o = sqrt(z_km_ω_5_to_1_1/y_km_1);
SIL = V_nom ^2/Z_o;
p_load = real(SIL)/100;
q_load = imag(SIL)/100;

l_seg = 10.0; # km 

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 0.1,
)

dommel_M1 = ETL.ExpParams(
    nothing, # N 
    1, # M
    100, # base line length 
    l_seg, 
    Z_c_abs_1, 
    z_km_1,
    y_km_1,
    z_km_ω_1,
    z_km_ω_5_to_1_1,
    Z_c_5_to_1_abs_1,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    p_load,
    q_load,
    1.0, # line_scale
    1.0 # load_scale
);

dommel_M3 = ETL.ExpParams(
    nothing, # N
    3, # M
    100, # base line length 
    l_seg, 
    Z_c_abs_3, 
    z_km_3,
    y_km_3,
    z_km_ω_3,
    z_km_ω_5_to_1_3,
    Z_c_5_to_1_abs_3,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    p_load,
    q_load,
    1.0, # line_scale
    1.0 # load_scale
);



ognjen_M1 = ETL.ExpParams(
    nothing, # N 
    1, # M
    100, # base line length 
    l_seg, 
    Z_c_abs_1, 
    z_km_1,
    y_km_1,
    z_km_ω_1,
    z_km_ω_5_to_1_1,
    Z_c_5_to_1_abs_1,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    p_load,
    q_load,
    1.0, # line_scale
    1.0 # load_scale
);
# Update
ognjen_M1.z_km_ω_5_to_1 = 0.01 + im*0.1;
ognjen_M1.z_km_ω = 0.01 + im*0.1;
ognjen_M1.y_km = im*10e-9*2*pi*60;
ognjen_M1.Z_c_5_to_1_abs = abs(sqrt(ognjen_M1.z_km_ω_5_to_1/ognjen_M1.y_km));
ognjen_M1.Z_c_abs = ognjen_M1.Z_c_5_to_1_abs;
ognjen_M1.line_scale = 0.5; # 50 km
ognjen_M1.l_seg = 5; 
ognjen_M1.p_load = 1.0;
ognjen_M1.q_load = 0.05;


function build_sim_from_file(file_name::String, dyn_lines::Bool, multi_segment::Bool, p::ExpParams, load_bus::String, inv_share, rating_scale)

    # build system
    sys = System(joinpath(pwd(), file_name));

    # Simulation time span
    tspan = (0.0, p.sim_params.t_max);
    perturbation = ETL.choose_disturbance(sys, p.perturbation, p);

    # Add load
    load = StandardLoad(
        name = "load1",
        available = true,
        bus = get_component(Bus, sys, load_bus),
        base_power = 100.0,
        constant_active_power = 0.0,
        constant_reactive_power = 0.0,
        impedance_active_power = p.p_load*p.load_scale,
        impedance_reactive_power = p.q_load*p.load_scale,
        current_active_power = 0.0,
        current_reactive_power = 0.0,
        max_constant_active_power = 0.0,
        max_constant_reactive_power = 0.0,
        max_impedance_active_power = p.p_load*p.load_scale,
        max_impedance_reactive_power = p.q_load*p.load_scale,
        max_current_active_power = 0.0,
        max_current_reactive_power = 0.0,
    );
    add_component!(sys, load);

    for g in get_components(PSY.Generator, sys)
        set_rating!(g, g.rating * rating_scale) # adjust the ratings of all generators
        set_active_power!(g, g.rating * inv_share); # adjust the share of all generators 
        # if g.bus.name == inv_bus # at the inverter bus 
        #     set_active_power!(g, g.rating * inv_share); # adjust P setpoint to be correct fraction of rating 
        # end
        
    end
    
    # build segments model
    if (multi_segment == true)
        sys = ETL.build_seg_model!(sys, p, dyn_lines, "");
    else
        sys = ETL.build_new_impedance_model!(sys, p, dyn_lines, "");
    end

    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);
    return sim 

end