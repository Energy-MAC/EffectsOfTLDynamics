module EffectsOfTLDynamics

using PowerSystems
using PowerSimulationsDynamics
using Sundials
using Plots
using PowerNetworkMatrices
using SparseArrays
using OrdinaryDiffEq


const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

function choose_disturbance(sys, dist)
    if dist == "CRC"
        # Control reference change
        g = get_component(DynamicInverter, sys, "generator-103-1")
        crc = ControlReferenceChange(0.25, g, :V_ref, 0.95)
    elseif dist == "NetworkSwitch"
        # NetworkSwitch
        yb = Ybus(sys).data
        new_yb = yb/1.0
        ns = NetworkSwitch(0.25, new_yb)
    elseif dist == "InfBusChange"
        # Source bus voltage change
        s_device = first(get_components(Source, sys))
        s_change = SourceBusVoltageChange(0.25, s_device, :V_ref, 1.048)
    else
        return error("Unknown disturbance")
    end
end

function build_sim(sys, tspan, perturbation, dyn_lines)
    sim = PSID.Simulation(
           MassMatrixModel, #Type of model used
           sys, #system
           pwd(), #folder to output results
           tspan, #time span
           perturbation, #Type of perturbation
           all_lines_dynamic = dyn_lines
       )
    return sim
end

function execute_sim(sim, p)
    exec = PSID.execute!(
           sim, #simulation structure
           Rodas4(), #Sundials DAE Solver
           dtmax = 1e-4, #0.02, #Arguments: Maximum timestep allowed
           abstol = p.abstol,
           maxiters = p.maxiters,
       );
       return exec
end

function results_sim(sim)
    results = read_results(sim)
    return results
end

function build_new_impedance_model(sys, p::ExpParams)
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
    y_ll = y_km_pu/2*l*(tanh(γ*l)/(γ*l))

        for l in get_components(Line, sys)
            l.r = real(z_ll)
            l.x = imag(z_ll)
            l.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
        end
    return sys
end

function build_seg_model(sys_segs, p::ExpParams)
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
    y_ll = y_km_pu/2*l*(tanh(γ*l)/(γ*l))
    N = p.N
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
    
    return sys_segs
end

function run_experiment(file_name, t_max, dist, line_model, p::ExpParams)
    # build system
    sys = System(joinpath(pwd(), file_name));

    # Simulation time span
    tspan = (0.0, t_max)

    # "CRC"
    # "NetworkSwitch"
    # "InfBusChange"
    perturbation = choose_disturbance(sys, dist)

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
        sys = build_seg_model(sys, p)
    else
        sys = build_new_impedance_model(sys, p)
    end
    # build simulation
    sim = build_sim(sys, tspan, perturbation, dyn_lines)
    show_states_initial_value(sim)
    # execute simulation
    exec = execute_sim(sim, p)
    # read results
    results = results_sim(sim)
    return results, sim
end

end