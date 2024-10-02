using PowerSystems
using PowerSimulationsDynamics
using Sundials
using PowerNetworkMatrices
using SparseArrays
using OrdinaryDiffEq
using CSV
using DataFrames

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
    elseif dist == "BranchTrip"
        dist_struct = BranchTrip(pp.t_fault, Line, pp.branch_trip_params.line_to_trip)
    else
        return error("Unknown disturbance")
    end
    return dist_struct
end

function build_sim(sys, tspan::Tuple{Float64, Float64}, perturbation, dyn_lines::Bool, p::ExpParams)
    if p.sim_params.solver == "Rodas4" || p.sim_params.solver == "FBDF"
        model = MassMatrixModel
    elseif p.sim_params.solver == "IDA"
        model = ResidualModel
    else
        return error("Unknown solver")
    end
    show_components(sys, Line)
    sim = PSID.Simulation(
           model, #Type of model used
           sys, #system
           pwd(), #folder to output results
           tspan, #time span
           perturbation #, #Type of perturbation
           # all_lines_dynamic = dyn_lines
       )
    return sim
end

function execute_sim!(sim, p::ExpParams)
    if p.sim_params.solver == "Rodas4"
        solver = Rodas4()
    elseif p.sim_params.solver == "IDA"
        solver = IDA()
    elseif p.sim_params.solver == "FBDF"
        solver = FBDF()
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

function build_new_impedance_model!(sys, p::ExpParams, dyn_lines::Bool, alg_line_name::String)

    Z_c_abs = p.Z_c_abs # Ω
    z_km = p.z_km # Ω/km
    y_km = p.y_km # S/km
    z_km_ω = p.z_km_ω # Ω/km
    z_km_ω_5_to_1 = p.z_km_ω_5_to_1
    Z_c_5_to_1_abs = p.Z_c_5_to_1_abs
    
    # z_km_pu = z_km/Z_c_abs
    z_km_ω_pu = z_km_ω/Z_c_abs
    y_km_pu = y_km*Z_c_abs
    γ = sqrt(z_km_ω*y_km)
    
    M = p.M
    if M == 1
        z_km_ω_pu = z_km_ω_5_to_1/Z_c_5_to_1_abs
        y_km_pu = y_km*Z_c_5_to_1_abs
        γ = sqrt(z_km_ω_5_to_1*y_km)
    end

    line_scale = p.line_scale

    for ll in get_components(Line, sys)
        # l = p.l
        l = p.l_dict[ll.name]*line_scale #km
        # println(l)
        z_ll = z_km_ω_pu*l*(sinh(γ*l)/(γ*l))
        y_ll = y_km_pu*l*(tanh(γ*l/2)/(γ*l/2))
        # println(z_ll)
        # println(y_ll)
        ll.r = real(z_ll)
        ll.x = imag(z_ll)
        ll.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
        if (ll.name != alg_line_name) && dyn_lines
            dyn_branch = DynamicBranch(get_component(Line, sys, ll.name))
            add_component!(sys, dyn_branch)
        end
    end

    return sys
end

function build_seg_model!(sys_segs, p::ExpParams, dyn_lines::Bool, alg_line_name::String)
    Z_c_abs = p.Z_c_abs # Ω
    z_km = p.z_km # Ω/km
    y_km = p.y_km # S/km
    z_km_ω = p.z_km_ω # Ω/km

    z_km_ω_5_to_1 = p.z_km_ω_5_to_1
    Z_c_5_to_1_abs = p.Z_c_5_to_1_abs
    
    z_km_pu = z_km/Z_c_abs
    y_km_pu = y_km*Z_c_abs
    z_km_ω_pu = z_km_ω/Z_c_abs
    γ = sqrt(z_km_ω*y_km)

    # N = p.N
    M = p.M
    l_seg = p.l_seg

    line_scale = p.line_scale

    if M == 1
        z_km_pu = z_km_ω_5_to_1/Z_c_5_to_1_abs
        y_km_pu = y_km*Z_c_5_to_1_abs
        z_km_ω_pu = z_km_ω_5_to_1/Z_c_5_to_1_abs
    end
    
    for ll in collect(get_components(Line, sys_segs))
        if ll.name == alg_line_name
            # l = p.l
            l = p.l_dict[ll.name]*line_scale #km
            # println(l)
            z_ll = z_km_ω_pu*l*(sinh(γ*l)/(γ*l))
            y_ll = y_km_pu*l*(tanh(γ*l/2)/(γ*l/2))
            # println(z_ll)
            # println(y_ll)
            # error("dayumn")
            ll.r = real(z_ll)
            ll.x = imag(z_ll)
            ll.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
            continue
        end
        l = p.l_dict[ll.name]*line_scale #km

        N = Int(ceil(l/l_seg))
        # l_seg = l/N
        z_seg_pu = z_km_pu*l_seg
        y_seg_pu = y_km_pu*l_seg
        bus_from = ll.arc.from
        bus_to = ll.arc.to
        # Create a bunch of Bus
        start_bus = bus_from
        for b_ix in 1 : N - 1
            println(b_ix)
            bus_to_create = ACBus(
                number = 1000000000 + 100000*bus_from.number + 100*bus_to.number + b_ix,
                name = bus_from.name * "-" * bus_to.name * "-internal-bus_" * string(b_ix),
                bustype = ACBusTypes.PQ,
                angle = bus_from.angle,
                magnitude = bus_from.magnitude,
                voltage_limits = bus_from.voltage_limits,
                base_voltage = bus_from.base_voltage,
                area = bus_from.area,
                load_zone = bus_from.load_zone
            )
            add_component!(sys_segs, bus_to_create)
            end_bus = bus_to_create
            for m in 1 : M
                line_to_create = Line(
                    name = ll.name * "_segment_" * string(b_ix) * "_branch_" * string(m),
                    available = true,
                    active_power_flow = ll.active_power_flow,
                    reactive_power_flow = ll.reactive_power_flow,
                    arc = Arc(from = start_bus, to = end_bus),
                    r = real(z_seg_pu[m]),
                    x = imag(z_seg_pu[m]),
                    b = (from = imag(y_seg_pu)/(2*M), to = imag(y_seg_pu)/(2*M)),
                    rate = ll.rate,
                    angle_limits = ll.angle_limits,
                )
                add_component!(sys_segs, line_to_create)
                if dyn_lines
                    dyn_branch = DynamicBranch(get_component(Line, sys_segs, line_to_create.name))
                    add_component!(sys_segs, dyn_branch)
                end
            end
            # println(get_name(line_to_create))
            # println("Bus From: $(line_to_create.arc.from.name), Bus To: $(line_to_create.arc.to.name)")
            start_bus = end_bus
        end
        for m in 1 : M
            line_to_create = Line(
                    name = ll.name * "_segment_" * string(N) * "_branch_" * string(m),
                    available = true,
                    active_power_flow = ll.active_power_flow,
                    reactive_power_flow = ll.reactive_power_flow,
                    arc = Arc(from = start_bus, to = bus_to),
                    r = real(z_seg_pu[m]),
                    x = imag(z_seg_pu[m]),
                    b = (from = imag(y_seg_pu)/(2*M), to = imag(y_seg_pu)/(2*M)),
                    rate = ll.rate,
                    angle_limits = ll.angle_limits,
            )
            add_component!(sys_segs, line_to_create)
            if dyn_lines
                dyn_branch = DynamicBranch(get_component(Line, sys_segs, line_to_create.name))
                add_component!(sys_segs, dyn_branch)
            end
        end
        # println(get_name(line_to_create))
        # println("Bus From: $(line_to_create.arc.from.name), Bus To: $(line_to_create.arc.to.name)")
        remove_component!(sys_segs, ll)
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

    if length(get_components(Bus, sys)) == 2
        ll = first(get_components(Line, sys))
        ll_alg = Line(
                    name = ll.name * "_static",
                    available = true,
                    active_power_flow = ll.active_power_flow,
                    reactive_power_flow = ll.reactive_power_flow,
                    arc = ll.arc,
                    r = ll.r,
                    x = ll.x,
                    b = (from = ll.b.from, to = ll.b.to),
                    rate = ll.rate,
                    angle_limits = ll.angle_limits,
            )
        # device = first(get_components(Generator, sys))
        # set_active_power!(device, 0.8)
        alg_line_name = ll_alg.name
        add_component!(sys, ll_alg)
        load = StandardLoad(
            name = "load1",
            available = true,
            bus = get_component(Bus, sys, "BUS 2"),
            base_power = 100.0,
            constant_active_power = 0.0,
            constant_reactive_power = 0.0,
            impedance_active_power = p.p_load,
            impedance_reactive_power = p.q_load,
            current_active_power = 0.0,
            current_reactive_power = 0.0,
            max_constant_active_power = 0.0,
            max_constant_reactive_power = 0.0,
            max_impedance_active_power = p.p_load,
            max_impedance_reactive_power = p.q_load,
            max_current_active_power = 0.0,
            max_current_reactive_power = 0.0,
        )
        add_component!(sys, load)
    end
    
    if length(get_components(Bus, sys)) == 9
        alg_line_name = p.perturbation_params.branch_trip_params.line_to_trip
    end

    for l in get_components(PSY.StandardLoad, sys)
        transform_load_to_constant_impedance(l)
        l.impedance_active_power = l.impedance_active_power * p.load_scale 
        l.impedance_reactive_power = l.impedance_reactive_power * p.load_scale 
    end
    
     for g in get_components(PSY.Generator, sys)
        if g.bus.bustype == ACBusTypes.PV
            # set_base_power!(g, g.base_power * p.load_scale)
            set_active_power!(g, g.active_power * p.load_scale)
            set_reactive_power!(g, g.reactive_power * p.load_scale)
            # set_rating!(g, g.rating * p.load_scale)
            # set_active_power_limits!(g, (min = g.active_power_limits.min * p.load_scale, max = g.active_power_limits.max * p.load_scale ))
            # set_reactive_power_limits!(g, (min = g.reactive_power_limits.min * p.load_scale, max = g.reactive_power_limits.max * p.load_scale ))
            # set_ramp_limits!(g, (up = g.ramp_limits.up * p.load_scale, down = g.ramp_limits.down * p.load_scale))      
        end
     end
    
    # build segments model
    if (multi_segment == true)
        sys = build_seg_model!(sys, p, dyn_lines, alg_line_name)
    else
        sys = build_new_impedance_model!(sys, p, dyn_lines, alg_line_name)
    end
    # build simulation

    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);
    show_states_initial_value(sim)

    # execute simulation
    exec = execute_sim!(sim, p);
    
    # read results
    results = results_sim(sim);
    return results, sim
end

function get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)
    df_imp = CSV.read(impedance_csv, DataFrame);
    c_km = CSV.read(capacitance_csv, DataFrame)."C"

    r_km = vec(zeros(M, 1))
    l_km = vec(zeros(M, 1))
    
    for m in 1 : M
        r_km[m,1] = df_imp[M, "R"*string(m)]
        l_km[m,1] = df_imp[M, "L"*string(m)]
    end
    
    f = 60
    ω = 2*pi*f

    x_km = ω*l_km
    z_km = (r_km + im*x_km)*factor_z
    g_km = [0.0]
    b_km = ω*c_km
    y_km = (g_km + im*b_km)*factor_y
    
    Y_ = 0
    for i in 1:M
        Y_ = Y_ + 1/(z_km[i])
    end

    z_km_ω = 1/Y_
    Z_c = sqrt(z_km_ω/y_km[1])

    r_km_5_to_1 = vec(zeros(5, 1))
    l_km_5_to_1 = vec(zeros(5, 1))
    T = 3
    for m in 1 : T
        r_km_5_to_1[m,1] = df_imp[T, "R"*string(m)]
        l_km_5_to_1[m,1] = df_imp[T, "L"*string(m)]
    end
    
    x_km_5_to_1 = ω*l_km_5_to_1
    z_km_5_to_1 = (r_km_5_to_1 + im*x_km_5_to_1)*factor_z

    Y_5_to_1 = 0 
    for i in 1:T
        Y_5_to_1 = Y_5_to_1 + 1/(z_km_5_to_1[i])
    end

    z_km_ω_5_to_1 = 1/Y_5_to_1
    Z_c_5_to_1 = sqrt(z_km_ω_5_to_1/y_km[1])
    
    return z_km, y_km[1], abs(Z_c), z_km_ω, z_km_ω_5_to_1, abs(Z_c_5_to_1)
end

# Verifying
function verifying(file_name, M, impedance_csv, capacitance_csv, p, factor_z, factor_y)
    sys = System(joinpath(pwd(), file_name))
    
    for ll in get_components(Line, sys)
        println("\n"*ll.name)
        println("z_pu_ll = $(ll.r) + j $(ll.x)")
        println("b_pu_ll = $(ll.b.from + ll.b.to)")

        for m in 1:M
            z_km, y_km, Z_c_abs, z_km_ω = get_line_parameters(impedance_csv, capacitance_csv, m, factor_z, factor_y)
            z_km_pu = z_km/Z_c_abs
            y_km_pu = y_km*Z_c_abs
            z_km_ω_pu = z_km_ω/Z_c_abs
            
            l = p.l_dict[ll.name]*p.line_scale
            println("Parallel branch impedances, with N = 1 and M = $(m)")
            z_pu_ll = z_km_pu*l
            #println("z_pu_ll = $(z_pu_ll)")
    
            b_pu_ll = imag(y_km_pu*l)
            println("b_pu_ll = $(b_pu_ll)")
            
            # z_pu_w_ll = z_km_ω_pu * l
            # println("z_pu_w_ll = $(z_pu_w_ll)")

            println("z_km_pu_ll = $(z_km_pu*l)")
    
            # abs_z_pu_w_ll = abs(z_pu_w_ll)
            # println("z_pu_w_ll = $(abs_z_pu_w_ll)")
        end
    end
end


function build_2bus_sim_from_file(file_name::String, dyn_lines::Bool, multi_segment::Bool, p::ExpParams, load_bus::String)
    # build system
    sys = System(joinpath(pwd(), file_name));

    # Simulation time span
    tspan = (0.0, p.sim_params.t_max);
    perturbation = choose_disturbance(sys, p.perturbation, p);
    
    # add extra algebraic line 
    
    ll = first(get_components(Line, sys))
    
    ll_alg = Line(
                name = ll.name * "_static",
                available = true,
                active_power_flow = ll.active_power_flow,
                reactive_power_flow = ll.reactive_power_flow,
                arc = ll.arc,
                r = ll.r,
                x = ll.x,
                b = (from = ll.b.from, to = ll.b.to),
                rate = ll.rate,
                angle_limits = ll.angle_limits,
        )
    alg_line_name = ll_alg.name
    add_component!(sys, ll_alg)

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

    # scale generation

    for g in get_components(PSY.Generator, sys)
        if g.bus.bustype == BusTypes.PV
            set_active_power!(g, g.active_power * p.load_scale);
            set_reactive_power!(g, g.reactive_power * p.load_scale);
        end
    end
    

    # build segments model
    if (multi_segment == true)
        sys = build_seg_model!(sys, p, dyn_lines, alg_line_name);
    else
        sys = build_new_impedance_model!(sys, p, dyn_lines, alg_line_name);
    end

    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);
    return sim 
end


function build_9bus_sim_from_file(file_name::String, dyn_lines::Bool, multi_segment::Bool, p::ExpParams)

    # build system
    sys = System(joinpath(pwd(), file_name));

    # Simulation time span
    tspan = (0.0, p.sim_params.t_max)
    perturbation = choose_disturbance(sys, p.perturbation, p)

    alg_line_name = p.perturbation_params.branch_trip_params.line_to_trip
    
    if (multi_segment == true)
        sys = build_seg_model!(sys, p, dyn_lines, alg_line_name)
    else
        sys = build_new_impedance_model!(sys, p, dyn_lines, alg_line_name)
    end

    # Scale loads according to load scale 
    for l in get_components(PSY.StandardLoad, sys)
        transform_load_to_constant_impedance(l)
        l.impedance_active_power = l.impedance_active_power * p.load_scale 
        l.impedance_reactive_power = l.impedance_reactive_power * p.load_scale 
    end

    # Scale generation by the load scale 
    for g in get_components(PSY.Generator, sys)
        set_active_power!(g, g.active_power * p.load_scale)
        set_reactive_power!(g, g.reactive_power * p.load_scale)
    end

    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);
    return sim 

end

function store_bus_voltages(res::PSID.SimulationResults, system::System, path::String)
    bus_numbers = sort(get_number.(get_components(Bus, system)))
    time, _ = get_voltage_magnitude_series(res, first(bus_numbers))
    df_voltages = DataFrame()
    # Store Time
    df_voltages[!, "Time"] = time
    for bus_n in bus_numbers
        _, bus_mag = get_voltage_magnitude_series(res, bus_n)
        _, bus_ang = get_voltage_angle_series(res, bus_n)
        df_voltages[!, "Bus-$(bus_n)-Mag"] = bus_mag
        df_voltages[!, "Bus-$(bus_n)-Ang"] = bus_ang
    end
    CSV.write(path, df_voltages, force = true)
end

function store_filter_currents(res::PSID.SimulationResults, sys::System, path::String)
    time, _ = get_state_series(res, (get_name(first(get_components(DynamicInverter, sys))), :ir_cnv))
    df_currents = DataFrame()
    # Store Time
    df_currents[!, "Time"] = time
    invs = collect(get_components(DynamicInverter, sys))
    
    for i in 1:length(invs)
        inv_name = invs[i].name
        _, ir_cnv = get_state_series(res, (inv_name, :ir_cnv))
        _, ii_cnv = get_state_series(res, (inv_name, :ii_cnv))
        df_currents[!, "$(inv_name)_ir_cnv"] = ir_cnv
        df_currents[!, "$(inv_name)_ii_cnv"] = ii_cnv
    end

    function store_tx_currents(res::PSID.SimulationResults, sys::System, path::String)
        time, _ = get_state_series(res, (get_name(first(get_components(DynamicInverter, sys))), :ir))
        df_currents = DataFrame()
        # Store Time
        df_currents[!, "Time"] = time
        invs = collect(get_components(DynamicInverter, sys))
        
        for i in 1:length(invs)
            inv_name = invs[i].name
            _, ir = get_state_series(res, (inv_name, :ir))
            _, ii = get_state_series(res, (inv_name, :ii))
            df_currents[!, "$(inv_name)_ir"] = ir
            df_currents[!, "$(inv_name)_ii"] = ii
            df_currents[!, "$(inv_name)_i_mag"] = sqrt(ir.^2+ ii.^2)
        end
    
        CSV.write(path, df_currents, force = true)
    end

    CSV.write(path, df_currents, force = true)
end

function store_stator_currents(res::PSID.SimulationResults, sys::System, path::String)
    time, _ = get_state_series(res, (get_name(first(get_components(DynamicGenerator, sys))), :ψ_d))
    df_currents = DataFrame()
    # Store Time
    df_currents[!, "Time"] = time
    gens = collect(get_components(DynamicGenerator, sys))
    
    for i in 1:length(gens)
        gen_name = gens[i].name
        _, ψ_d = get_state_series(res, (gen_name, :ψ_d))
        _, ψ_q = get_state_series(res, (gen_name, :ψ_q))
        _, ed_pp = get_state_series(res, (gen_name, :ed_pp))
        _, eq_pp = get_state_series(res, (gen_name, :eq_pp))

        xd_pp = get_Xd_pp(gens[i])
        xq_pp = get_Xq_pp(gens[i])
        
        id = 1/xd_pp * (eq_pp - ψ_d)
        iq = 1/xq_pp * (- ed_pp - ψ_q)
        
        df_currents[!, "$(gen_name)_id"] = id
        df_currents[!, "$(gen_name)_iq"] = iq
        df_currents[!, "$(gen_name)_i_mag"] = sqrt(id.^2+ iq.^2)
    end

    CSV.write(path, df_currents, force = true)
end

function store_branch_power_flows(res::PSID.SimulationResults, sys::System, path::String)
    time, _ = get_activepower_branch_flow(res, get_name(first(get_components(Line, sys))), :from)
    df_power = DataFrame()
    # Store Time
    df_power[!, "Time"] = time
    lines = collect(get_components(Line, sys))
    
    for i in 1:length(lines)
        line_name = lines[i].name
        _, P = get_activepower_branch_flow(res, line_name, :from)
        _, Q = get_reactivepower_branch_flow(res, line_name, :from)
        df_power[!, "$(line_name)_P"] = P
        df_power[!, "$(line_name)_Q"] = Q
    end

    CSV.write(path, df_power, force = true)
end

function store_generator_speeds(res::PSID.SimulationResults, sys::System, path::String)
    time, _ = get_state_series(res, (get_name(first(get_components(DynamicGenerator, sys))), :ω))
    df_speeds = DataFrame()
    # Store Time
    df_speeds[!, "Time"] = time
    gens = collect(get_components(DynamicGenerator, sys))
    
    for i in 1:length(gens)
        gen_name = gens[i].name
        _, w = get_state_series(res, (gen_name, :ω))
        df_speeds[!, "$(gen_name)_ω"] = w
    end
    CSV.write(path, df_speeds)
end

function store_branch_currents(res::PSID.SimulationResults, sys::System, path::String)
    time, _ = get_real_current_branch_flow(res, get_name(first(get_components(Line, sys))))
    df_currents = DataFrame()
    # Store Time
    df_currents[!, "Time"] = time
    lines = collect(get_components(Line, sys))
    
    for i in 1:length(lines)
        line_name = lines[i].name
        _, I_R = get_real_current_branch_flow(res, line_name)
        _, I_I = get_imaginary_current_branch_flow(res, line_name)
        df_currents[!, "$(line_name)_I_R"] = I_R
        df_currents[!, "$(line_name)_I_I"] = I_I
    end

    CSV.write(path, df_currents, force = true)
end

function extract_simulation_data(sim, dt=nothing)
    # Initialize the dictionary to hold all the data
    res_dict = Dict()
    # Iterate over each key in the global index (assuming these are your device names like "V_1", "generator-102-1")
    for (device_name, vars) in sim.results.global_index
        # Ensure each device has its own sub-dictionary
        res_dict[device_name] = Dict()
        # Iterate over each variable and its index within the device
        for (var_name, index) in vars
            # Get the time series for the specific variable at all times
            ts, state = PSID.get_state_series(sim.results, (device_name, var_name), dt=dt)  # Assuming you have a function that can do this
            # Assign the time series to the correct variable in the nested dictionary
            res_dict[device_name][var_name] = (ts, state)
        end
    end

    return res_dict
end

export choose_disturbance
export build_sim
export execute_sim
export results_sim
export build_new_impedance_model
export build_seg_model
export run_experiment
export get_line_parameters
export verifying
export store_bus_voltages
export store_filter_currents
export store_branch_power_flows
export store_generator_speeds
export store_branch_currents
export extract_simulation_data
