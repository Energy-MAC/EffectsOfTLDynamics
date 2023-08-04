using PowerSystems
using PowerSimulationsDynamics
using Sundials
using PowerNetworkMatrices
using SparseArrays
using OrdinaryDiffEq

include("ExperimentStructs.jl")
using .ExperimentStructs

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

function choose_disturbance(sys, dist::String, p::ExpParams)
    pp = p.perturbation_params
    
    if dist == "BIC"
        dist_struct = BranchImpedanceChange(pp.t_fault, Line, pp.branch_impedance_change_params.line_name, pp.branch_impedance_change_params.multiplier)
    elseif dist == "GenTrip"
        g = get_component(pp.gen_trip_params.source_type, sys, pp.gen_trip_params.gen_name)
        dist_struct = GeneratorTrip(pp.t_fault, g)
    elseif dist == "CRC"
        # Control reference change
        g = get_component(pp.crc_params.source_type, sys, pp.crc_params.gen_name)
        dist_struct = ControlReferenceChange(pp.t_fault, g, pp.crc_params.var_to_change, pp.crc_params.ref_value)
    elseif dist == "LoadChange"
        # Load change
        l = get_component(pp.load_change_params.load_type, sys, pp.load_change_params.load_name)
        dist_struct = LoadChange(pp.t_fault, l, pp.load_change_params.var_to_change, pp.load_change_params.ref_value)
    elseif dist == "LoadTrip"
        # Load trip
        l = get_component(pp.load_trip_params.load_type, sys, pp.load_trip_params.load_name)
        dist_struct = LoadTrip(pp.t_fault, l)
    elseif dist == "InfBusChange"
        # Source bus voltage change
        s_device = first(get_components(Source, sys))
        dist_struct = SourceBusVoltageChange(pp.t_fault, s_device, pp.source_bus_voltage_change_params.var_to_change, pp.source_bus_voltage_change_params.ref_value)
    else
        return error("Unknown disturbance")
    end
    return dist_struct
end

function build_sim(sys, tspan::Tuple{Float64, Float64}, perturbation, dyn_lines::Bool, p::ExpParams)
    if p.sim_params.solver == "Rodas4"
        model = MassMatrixModel
    elseif p.sim_params.solver == "IDA"
        model = ResidualModel
    else
        return error("Unknown solver")
    end
    sim = PSID.Simulation(
           model, #Type of model used
           sys, #system
           pwd(), #folder to output results
           tspan, #time span
           perturbation, #Type of perturbation
           all_lines_dynamic = dyn_lines
       )
    return sim
end

function execute_sim!(sim, p::ExpParams)
    if p.sim_params.solver == "Rodas4"
        solver = Rodas4()
    elseif p.sim_params.solver == "IDA"
        solver = IDA()
    else
        return error("Unknown solver")
    end
    exec = PSID.execute!(
           sim, #simulation structure
           solver, #Sundials DAE Solver
           dtmax = p.sim_params.dtmax, #0.02, #Arguments: Maximum timestep allowed
           abstol = p.sim_params.abstol,
           maxiters = p.sim_params.maxiters,
       );
       return exec
end

function results_sim(sim)
    results = read_results(sim)
    return results
end

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
    y_ll = y_km_pu/2*l*(tanh(γ*l/2)/(γ*l/2))

        for l in get_components(Line, sys)
            l.r = real(z_ll)
            l.x = imag(z_ll)
            l.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
        end
    return sys
end

function build_seg_model!(sys_segs, p::ExpParams)
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
    y_ll = y_km_pu/2*l*(tanh(γ*l/2)/(γ*l/2))
    N = p.N
    l_prime = l/N

    z_seg_pu = z_km_pu*l_prime*(sinh(γ*l_prime)/(γ*l_prime))
    y_seg_pu = y_km_pu/2*l_prime*(tanh(γ*l_prime/2)/(γ*l_prime/2))

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

function run_experiment(file_name::String, line_model::String, p::ExpParams)
    # build system
    sys = System(joinpath(pwd(), file_name));

    # Simulation time span
    tspan = (0.0, p.sim_params.t_max)

    # "CRC"
    # "NetworkSwitch"
    # "InfBusChange"
    perturbation = choose_disturbance(sys, p.perturbation, p)

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
    show_states_initial_value(sim)
    # execute simulation
    exec = execute_sim!(sim, p)
    # read results
    results = results_sim(sim)
    return results, sim
end

export choose_disturbance
export build_sim
export execute_sim
export results_sim
export build_new_impedance_model
export build_seg_model
export run_experiment
