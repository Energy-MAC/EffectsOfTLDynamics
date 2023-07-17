using PowerSystems
using PowerSimulationsDynamics
using Sundials
using Plots
using PowerNetworkMatrices
using SparseArrays
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "test_sys.json"));
sys_dyn = deepcopy(sys)
sys_segs = deepcopy(sys)

# Simulation time span
tspan = (0.0, 30.0)

### Disturbances
# Function to choose disturbance

function choose_disturbance(dist)
    if dist == "CRC"
        # Control reference change
        g = get_component(DynamicInverter, sys, "generator-102-1")
        crc = ControlReferenceChange(1.0, g, :V_ref, 0.95)
    elseif dist == "NetworkSwitch"
        # NetworkSwitch
        yb = Ybus(sys).data
        new_yb = yb/1.0
        ns = NetworkSwitch(1.0, new_yb)
    elseif dist == "InfBusChange"
        # Source bus voltage change
        s_device = first(get_components(Source, sys))
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.048)
    else
        return error("Unknown disturbance")
    end
end

# "CRC"
# "NetworkSwitch"
# "InfBusChange"

dist = "CRC"
perturbation = choose_disturbance(dist)

# Build static simulation
sim = PSID.Simulation(
           ResidualModel, #Type of model used
           sys, #system
           pwd(), #folder to output results
           tspan, #time span
           perturbation, #Type of perturbation
       )
show_states_initial_value(sim)

PSID.execute!(
           sim, #simulation structure
           IDA(), #Sundials DAE Solver
           dtmax = 0.02, #Arguments: Maximum timestep allowed
       );

results = read_results(sim)

#####

# Build dynamic simulation
sim_dyn = PSID.Simulation(
           ResidualModel, #Type of model used
           sys_dyn, #system
           pwd(), #folder to output results
           tspan, #time span
           perturbation, #Type of perturbation
           all_lines_dynamic = true
       )
show_states_initial_value(sim_dyn)

PSID.execute!(
    sim_dyn, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02, #Arguments: Maximum timestep allowed
    abstol = 1.0e-5
);

results_dyn = read_results(sim_dyn)

#####
# Build dynamic segment simulation
sys_segs = deepcopy(sys)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
z_km = r_km + im*x_km # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
y_km = g_km + im*b_km
z_km_pu = z_km/Z_c
y_km_pu = y_km*Z_c
l = 100 #km
γ = sqrt(z_km*y_km)
z_ll = z_km_pu*l*(sinh(γ*l)/(γ*l))
y_ll = y_km_pu/2*l*(tanh(γ*l)/(γ*l))
N = 1
l_prime = l/N

z_seg_pu = z_km_pu*l_prime*(sinh(γ*l_prime)/(γ*l_prime))
y_seg_pu = y_km_pu/2*l_prime*(tanh(γ*l_prime)/(γ*l_prime))

for l in collect(get_components(Line, sys_segs))
    bus_from = l.arc.from
    bus_to = l.arc.to
    # Create a bunch of Bus
    start_bus = bus_from
    for b_ix in 1:N - 1
        println(b_ix)
        bus_to_create = Bus(
            number = 1000000000 + 100000*bus_from.number + 100*bus_to.number + b_ix,
            name = bus_from.name * "-" * bus_to.name * "-internal-bus_" * string(b_ix),
            bustype = BusTypes.PQ,
            angle = bus_from.angle,
            magnitude = bus_from.magnitude,
            voltage_limits = bus_from.voltage_limits,
            base_voltage = bus_from.base_voltage,
            area = bus_from.area,
            load_zone = bus_from.load_zone
        )
        add_component!(sys_segs, bus_to_create)
        end_bus = bus_to_create
        line_to_create = Line(
            name = l.name * "_segment_" * string(b_ix),
            available = true,
            active_power_flow = l.active_power_flow,
            reactive_power_flow = l.reactive_power_flow,
            arc = Arc(from = start_bus, to = end_bus),
            r = real(z_seg_pu),
            x = imag(z_seg_pu),
            b = (from = imag(y_seg_pu)/2, to = imag(y_seg_pu)/2),
            rate = l.rate,
            angle_limits = l.angle_limits,
        )
        add_component!(sys_segs, line_to_create)
        # println(get_name(line_to_create))
        # println("Bus From: $(line_to_create.arc.from.name), Bus To: $(line_to_create.arc.to.name)")
        start_bus = end_bus
    end
    line_to_create = Line(
            name = l.name * "_segment_" * string(N),
            available = true,
            active_power_flow = l.active_power_flow,
            reactive_power_flow = l.reactive_power_flow,
            arc = Arc(from = start_bus, to = bus_to),
            r = real(z_seg_pu),
            x = imag(z_seg_pu),
            b = (from = imag(y_seg_pu)/2, to = imag(y_seg_pu)/2),
            rate = l.rate,
            angle_limits = l.angle_limits,
    )
    add_component!(sys_segs, line_to_create)
    # println(get_name(line_to_create))
    # println("Bus From: $(line_to_create.arc.from.name), Bus To: $(line_to_create.arc.to.name)")
    remove_component!(sys_segs, l)
end

sim_segs = PSID.Simulation(
           ResidualModel, #Type of model used
           sys_segs, #system
           pwd(), #folder to output results
           tspan, #time span
           perturbation, #Type of perturbation
           all_lines_dynamic = true
       )
show_states_initial_value(sim_segs)

PSID.execute!(
    sim_segs, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02, #Arguments: Maximum timestep allowed
    abstol = 1.0e-8
);

results_segs = read_results(sim_segs)

# Plotting
vr = get_state_series(results, ("generator-102-1", :vr_filter));
plot(vr, xlabel = "time", ylabel = "vr p.u.", label = "vr")
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));
plot!(vr_dyn, xlabel = "time", ylabel = "vr", label = "vr_dyn")
vr_segs = get_state_series(results_segs, ("generator-102-1", :vr_filter));
plot!(vr_segs, xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(N)", xlims=(0.99, 1.01))
plot!(xlims=(0.99, 1.01))
p_branch = get_activepower_branch_flow(results, "BUS 1-BUS 2-i_1", :from)
p_branch_dyn = get_activepower_branch_flow(results_dyn, "BUS 1-BUS 2-i_1", :from)
p_branch_seg = get_activepower_branch_flow(results_segs, "BUS 1-BUS 2-i_1_segment_1", :from)

plot(p_branch)
plot!(p_branch_dyn)
plot!(p_branch_seg,xlims=(0.99,1.07))