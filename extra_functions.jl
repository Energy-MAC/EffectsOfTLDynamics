# Define function to build simulation using multi-seg dynamic lines

function get_ms_dyn_sim(file_name, t_max, dist, p::ExpParams)
    sys = System(joinpath(pwd(), file_name));
    tspan = (0.0, t_max)
    perturbation = choose_disturbance(sys, dist)
    sys = build_seg_model(sys, p)
    sim = build_sim(sys, tspan, perturbation, true);
#    show_states_initial_value(sim);
    return sim
end

function get_ms_dyn_sys(file_name, p::ExpParams)
    sys = System(joinpath(pwd(), file_name));
    sys = build_seg_model(sys, p)
    return sys
end


function plot_N_vs_second_λ_for_Ls(sys, Lrange, Nrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan);
    # Plots second largest eigenvalue as a function of N, looping through different line lengths. 
    plot()
    for l in Lrange;
        second_λ = [];
        for N in Nrange;
            p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
            og_sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
            sys_ms = build_seg_model(og_sys_copy, p)
            sim = build_sim(sys_ms, tspan, perturbation, true);
            ss = small_signal_analysis(sim)
            push!(second_λ, real(ss.eigenvalues[end-1]))
        end
        display(plot!(Nrange, second_λ, seriestype=:scatter, label=string(l)*"km"))
    end

    xlabel!("N segs")
    ylabel!("Second largest real λ")

end


function plot_all_eigs_for_Ns(sys, Nrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan);
    plot()

    for N in Nrange;
        p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
        og_sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
        sys_ms = build_seg_model(og_sys_copy, p)
        sim = build_sim(sys_ms, tspan, perturbation, true);
        ss = small_signal_analysis(sim)
        display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="N="*string(N)))
    end

    xlabel!("Real")
    ylabel!("Imag")

end

function make_sys_from_filename(file_name);
    sys = System(joinpath(pwd(), file_name));
    return sys
end

function plot_N_vs_xth_largest_λ_for_Ls(x, sys, Lrange, Nrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan);
    # Plots second largest eigenvalue as a function of N, looping through different line lengths. 
    plot()
    for l in Lrange;
        x_λ = [];
        for N in Nrange;
            p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
            og_sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
            sys_ms = build_seg_model(og_sys_copy, p)
            sim = build_sim(sys_ms, tspan, perturbation, true);
            ss = small_signal_analysis(sim)
            push!(x_λ, real(ss.eigenvalues[end-(x-1)]))
        end
        display(plot!(Nrange, x_λ, seriestype=:scatter, label=string(l)*"km"))
    end

    xlabel!("N segs")
    ylabel!("Largest non-zero real λ ("*string(x)*")")

end

function save_current_fig();
    savefig("/Users/ruthkravis/Downloads/plot.png")
end

function get_modified_twobus_sys();
    file_name = "test_sys.json";
    sys = System(joinpath(pwd(), file_name));

    inv = get_component(DynamicInverter, sys, "generator-102-1"); # Get first inverter
    inv2 = deepcopy(inv); # Create second inverter 
    inv2.name = "generator-101-1"; # Rename second inverter

    inv_static = get_component(StaticInjection, sys, "generator-102-1"); # Get static injection associated with first inverter
    inv2_static = deepcopy(inv_static); # Copy for second inverter
    inv2_static.name = "generator-101-1"; # Rename static injection to match dynamic component 
    inv2_static.dynamic_injector = inv2; # assign second inverter dynamic injector to second inverter static injection object 
    inv2_static.time_series_container = InfrastructureSystems.TimeSeriesContainer(); # To fix an error I was getting 

    inf_bus = get_component(Source, sys, "InfBus"); # Get the infinite bus 
    inv2_static.bus = inf_bus.bus; # Use the bus info from inf bus to update the bus of inv2 
    # inv2_static.bus.bustype = BusTypes.SLACK; # May need to change bus type?

    add_component!(sys, inv2_static); # add the second inverter 
    remove_component!(sys,inf_bus); # remove the infinite bus source 

    return sys 
end


function plot_nonzero_eig_fixed_Lseg_for_Ls(x, Lseg, sys, Lrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan);
    # Plots the largest real eigenvalue part that is not zero, given a fixed segment length Lseg, for different line lengths Lrange. 

    plot()
    x_λ = []
    for l in Lrange;
        N = l/Lseg;
        p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
        og_sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
        sys_ms = build_seg_model(og_sys_copy, p)
        sim = build_sim(sys_ms, tspan, perturbation, true);
        ss = small_signal_analysis(sim)
        push!(x_λ, real(ss.eigenvalues[end-(x-1)]))
        
    end
    plot!(Lrange, x_λ, seriestype=:scatter)

    xlabel!("Line length (km)")
    ylabel!("Largest non-zero real λ ("*string(x)*")")
    title!(string(Lseg)*" km segments")

end

