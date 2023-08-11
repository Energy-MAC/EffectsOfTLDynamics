using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

function build_new_impedance_model!(sys, p::ExpParams)
    Z_c = p.Z_c # Ω
    r_km = p.r_km # Ω/km
    x_km = p.x_km # Ω/km
    z_km = r_km + im*x_km # Ω/km
    
    g_km = p.g_km # S/km
    b_km = p.b_km # S/km
    y_km = g_km + im*b_km
    
    z_km_pu = z_km/Z_c
    y_km_pu = y_km*Z_c
    
    l = p.l #km
    γ = sqrt(z_km*y_km)
    z_ll = z_km_pu*l*(sinh(γ*l)/(γ*l))
    y_ll = y_km_pu*l*(tanh(γ*l/2)/(γ*l/2))

        for ll in get_components(Line, sys)
            ll.r = real(z_ll)
            ll.x = imag(z_ll)
            ll.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
        end
    return sys
end

sim_p = SimParams(
    abstol = 1e13,
    reltol = 1e10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 2.0,
)
# "OMIB.json"
# "9bus.json"
file_name = "OMIB.json"

# "BIC"
# "GenTrip"
# "CRC"
# "LoadChange"
# "LoadTrip"
# "InfBusChange"
perturbation = "CRC"

M = 1
# Z_c, r_km, x_km, g_km, b_km = get_line_parameters(data_file, M)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
l = 1000 #km
N = nothing
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation)
p = ExpParams(N, M, l, Z_c, r_km, x_km, g_km, b_km, sim_p, perturbation, perturbation_params)

sys = System(joinpath(pwd(), file_name));

# Simulation time span
tspan = (0.0, p.sim_params.t_max)

# "CRC"
# "NetworkSwitch"
# "InfBusChange"
perturbation = choose_disturbance(sys, p.perturbation, p)
line_model = "Algebraic"

# choose line model
if line_model == "Algebraic"
    # Algebraic Pi Lines
    dyn_lines = false
    multi_segment = false
elseif line_model == "Dynamic"
    # Dynamic Pi Lines
    dyn_lines = true
    multi_segment = false
elseif line_model == "Multi-Segment Algebraic"
    # Multi-Segment Algebraic Pi Lines
    dyn_lines = false
    multi_segment = true
elseif line_model == "Multi-Segment Dynamic"
    # Multi-Segment Dynamic Pi Lines
    dyn_lines = true
    multi_segment = true
else
    return error("Unknown line model")
end


# build segments model
if (multi_segment == true)
    sys = build_seg_model!(sys, p)
else
    sys = build_new_impedance_model!(sys, p)
end
# build simulation
sim = build_sim(sys, tspan, perturbation, dyn_lines, p)

s = small_signal_analysis(sim)