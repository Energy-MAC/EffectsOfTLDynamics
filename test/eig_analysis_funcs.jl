using PowerSystems
using PowerSimulationsDynamics
using Sundials
using PowerNetworkMatrices
using SparseArrays
using OrdinaryDiffEq
using CSV
using JLD2

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;


function save_max_nonzero_eig!(sim, output);
    if sim.status != PSID.BUILD_FAILED;
        ss = small_signal_analysis(sim);
        if maximum(real(ss.eigenvalues)) == 0.0;
            # Choose second largest eig 
            push!(output, real(ss.eigenvalues[end-1]));
        else
            push!(output, maximum(real(ss.eigenvalues)));
        end
    else
        push!(output, NaN)
    end
end


function generate_2bus_line_sweep_data(file_name::String, p1, p3, load_scale::Float64, line_scales::Vector, l_seg::Float64, load_bus::String)
    
    mssb_results = [];
    msmb_results = [];
    dyn_results = [];
    alg_results = [];

    p1.load_scale = load_scale;
    p3.load_scale = load_scale;

    p1.l_seg = l_seg;
    p3.l_seg = l_seg;

    for line_scale = line_scales;
        p1.line_scale = line_scale;
        p3.line_scale = line_scale;

        # Algebraic
        sim_alg = ETL.build_2bus_sim_from_file(file_name, false, false, p1, load_bus);
        save_max_nonzero_eig!(sim_alg, alg_results);

        # Dynamic pi  
        sim_dyn = ETL.build_2bus_sim_from_file(file_name, true, false, p1, load_bus);
        save_max_nonzero_eig!(sim_dyn, dyn_results);

        # MSSB
        sim_ms = ETL.build_2bus_sim_from_file(file_name, true, true, p1, load_bus);
        save_max_nonzero_eig!(sim_ms, mssb_results);

        # Do MSMB 
        sim_msmb = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus);
        save_max_nonzero_eig!(sim_msmb, msmb_results);

    end

    return alg_results, dyn_results, mssb_results, msmb_results

end

function generate_2bus_loading_v_boundary_data(file_name::String, p1, p3, load_scales::Vector, line_scales::Vector, l_seg::Float64, load_bus::String)

    alg_lims = [];
    dyn_lims = [];
    mssb_lims = [];
    msmb_lims = [];
    p1.l_seg = l_seg;
    p3.l_seg = l_seg;

    idx = 1;
    for load_scale in load_scales;
        mssb_results = [];
        msmb_results = [];
        dyn_results = [];
        alg_results = [];

        p1.load_scale = load_scale
        p3.load_scale = load_scale

        current_line_scales = line_scales[idx];
        for line_scale = current_line_scales;
            p1.line_scale = line_scale
            p3.line_scale = line_scale

            sim_alg = ETL.build_2bus_sim_from_file(file_name, false, false, p1, load_bus);
            save_max_nonzero_eig!(sim_alg, alg_results)

            # Dynamic pi  
            sim_dyn = ETL.build_2bus_sim_from_file(file_name, true, false, p1, load_bus);
            save_max_nonzero_eig!(sim_dyn, dyn_results)

            # MSSB
            sim_ms = ETL.build_2bus_sim_from_file(file_name, true, true, p1, load_bus);
            save_max_nonzero_eig!(sim_ms, mssb_results)

            # Do MSMB 
            sim_msmb = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus);
            save_max_nonzero_eig!(sim_msmb, msmb_results)

        end

        # Find first length where we lose stability - NaN indicates that we don't encounter stability limit 
        if findfirst(alg_results.>0) !== nothing;
            push!(alg_lims, current_line_scales[findfirst(alg_results.>0)])
        else
            push!(alg_lims, NaN)
        end
        if findfirst(dyn_results.>0) !== nothing;
            push!(dyn_lims, current_line_scales[findfirst(dyn_results.>0)])
        else
            push!(dyn_lims, NaN)
        end
        if findfirst(mssb_results.>0) !== nothing;
            push!(mssb_lims, current_line_scales[findfirst(mssb_results.>0)])
        else
            push!(mssb_lims, NaN)
        end
        if findfirst(msmb_results.>0) !== nothing;
            push!(msmb_lims, current_line_scales[findfirst(msmb_results.>0)])
        else
            push!(msmb_lims, NaN)
        end
        idx += 1
    end

    # plot1 = plot(load_scales, alg_lims,xlabel="Load scale", ylabel="Line length @ stability boundary", label="Algebraic", seriestype=:line, size=(800,600))
    # plot!(load_scales, dyn_lims, seriestype=:line, label="Dynpi")
    # plot!(load_scales, mssb_lims, seriestype=:line, linestyle=:dashdot, label="MSSB lseg max="*string(l_seg))
    # plot!(load_scales, msmb_lims, seriestype=:line, linestyle=:dash, label="MSMB lseg max="*string(l_seg))
    # title!("Base p="*string(round(p1.p_load, digits=2))*", q="*string(round(p1.q_load, digits=2)))

    return alg_lims, dyn_lims, mssb_lims, msmb_lims

end

function generate_9bus_line_sweep_data(file_name::String, p1, p3, load_scale::Float64, line_scales::Vector, l_seg::Float64)
    
    mssb_results = [];
    msmb_results = [];
    dyn_results = [];
    alg_results = [];

    p1.load_scale = load_scale;
    p3.load_scale = load_scale;

    p1.l_seg = l_seg;
    p3.l_seg = l_seg;

    for line_scale = line_scales;

        p1.line_scale = line_scale
        p3.line_scale = line_scale

        # Algebraic 
        sim_alg = ETL.build_9bus_sim_from_file(file_name, false, false, p1);
        save_max_nonzero_eig!(sim_alg, alg_results)

        # Dynamic pi  
        sim_dyn = ETL.build_9bus_sim_from_file(file_name, true, false, p1);
        save_max_nonzero_eig!(sim_dyn, dyn_results)

        # MSSB
        sim_ms = ETL.build_9bus_sim_from_file(file_name, true, true, p1);
        save_max_nonzero_eig!(sim_ms, mssb_results)

        # Do MSMB 
        sim_msmb = ETL.build_9bus_sim_from_file(file_name, true, true, p3);
        save_max_nonzero_eig!(sim_msmb, msmb_results)

    end

    return alg_results, dyn_results, mssb_results, msmb_results

end


function generate_9bus_individual_line_sweep_data(file_name::String, p1, p3, load_scale::Float64, line_scales::Vector, l_seg::Float64, line_to_scale::String)
    
    mssb_results = [];
    msmb_results = [];
    dyn_results = [];
    alg_results = [];

    p1.load_scale = load_scale;
    p3.load_scale = load_scale;

    p1.l_seg = l_seg;
    p3.l_seg = l_seg;

    original_line_length = p1.l_dict[line_to_scale];
    p1.line_scale = 1.0
    p3.line_scale = 1.0

    for line_scale = line_scales;

        # Choose single line to scale 
        p1.l_dict[line_to_scale] = original_line_length*line_scale
        p3.l_dict[line_to_scale] = original_line_length*line_scale

        # Algebraic 
        sim_alg = ETL.build_9bus_sim_from_file(file_name, false, false, p1);
        save_max_nonzero_eig!(sim_alg, alg_results)

        # Dynamic pi  
        sim_dyn = ETL.build_9bus_sim_from_file(file_name, true, false, p1);
        save_max_nonzero_eig!(sim_dyn, dyn_results)

        # MSSB
        sim_ms = ETL.build_9bus_sim_from_file(file_name, true, true, p1);
        save_max_nonzero_eig!(sim_ms, mssb_results)

        # Do MSMB 
        sim_msmb = ETL.build_9bus_sim_from_file(file_name, true, true, p3);
        save_max_nonzero_eig!(sim_msmb, msmb_results)

    end

    # Reset line length to original 
    p1.l_dict[line_to_scale] = original_line_length
    p3.l_dict[line_to_scale] = original_line_length

    return alg_results, dyn_results, mssb_results, msmb_results

end

function generate_9bus_load_scale_v_boundary_data(file_name::String, p1, p3, load_scales::Vector, line_scales::Vector, l_seg::Float64)

    alg_lims = [];
    dyn_lims = [];
    mssb_lims = [];
    msmb_lims = [];

    p1.l_seg = l_seg;
    p3.l_seg = l_seg;

    for load_scale in load_scales;
        mssb_results = [];
        msmb_results = [];
        dyn_results = [];
        alg_results = [];

        p1.load_scale = load_scale
        p3.load_scale = load_scale

        current_lscale = line_scales[idx][1];
        for line_scale = current_lscale;
    
            p1.line_scale = line_scale
            p3.line_scale = line_scale
    
            # Algebraic 
            sim_alg = ETL.build_9bus_sim_from_file(file_name, false, false, p1);
            save_max_nonzero_eig!(sim_alg, alg_results)
    
            # Dynamic pi  
            sim_dyn = ETL.build_9bus_sim_from_file(file_name, true, false, p1);
            save_max_nonzero_eig!(sim_dyn, dyn_results)
    
            # MSSB
            sim_ms = ETL.build_9bus_sim_from_file(file_name, true, true, p1);
            save_max_nonzero_eig!(sim_ms, mssb_results)
    
            # Do MSMB 
            sim_msmb = ETL.build_9bus_sim_from_file(file_name, true, true, p3);
            save_max_nonzero_eig!(sim_msmb, msmb_results)
    
        end
    
        # Find first length where we lose stability - NaN indicates that we don't encounter stability limit 
        if findfirst(alg_results.>0) !== nothing;
            push!(alg_lims, current_lscale[findfirst(alg_results.>0)])
        else
            push!(alg_lims, NaN)
        end
        if findfirst(dyn_results.>0) !== nothing;
            push!(dyn_lims, current_lscale[findfirst(dyn_results.>0)])
        else
            push!(dyn_lims, NaN)
        end
        if findfirst(mssb_results.>0) !== nothing;
            push!(mssb_lims, current_lscale[findfirst(mssb_results.>0)])
        else
            push!(mssb_lims, NaN)
        end
        if findfirst(msmb_results.>0) !== nothing;
            push!(msmb_lims, current_lscale[findfirst(msmb_results.>0)])
        else
            push!(msmb_lims, NaN)
        end
        idx += 1
    end
    
    return alg_lims, dyn_lims, mssb_lims, msmb_lims

end


function plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, p1, size::Tuple, linewidth::Int, ylims::Tuple, ylabel::String, xlabel::String, title::String, legend_on::Bool)

    plot1 = plot(load_scales, alg,xlabel=xlabel, ylabel=ylabel, label="Algebraic", seriestype=:line, linewidth=linewidth, size=size,legend=:outertopright)
    plot!(load_scales, dyn, seriestype=:line, linewidth=linewidth, label="Dynpi", legend=legend_on)
    plot!(load_scales, mssb, seriestype=:line, linewidth=linewidth, linestyle=:dashdot, label="MSSB lseg max="*string(l_seg))
    plot!(load_scales, msmb, seriestype=:line, linewidth=linewidth, linestyle=:dash, label="MSMB lseg max="*string(l_seg))
    title!(title)
    #title!("Base p="*string(round(p1.p_load, digits=2))*", q="*string(round(p1.q_load, digits=2)))
    ylims!(ylims)
    return plot1

end

function plot_2bus_line_sweep(alg, dyn, mssb, msmb, line_scales, l_seg, size, linewidth, title, legendloc)

    plot1 = plot(line_scales, alg,xlabel="Line scale (base=100km)", ylabel="Max real λ != 0", label="Algebraic", legend=legendloc, seriestype=:line, linewidth=linewidth, size=size)
    plot!(line_scales, dyn, seriestype=:line, linewidth=linewidth, label="Dynpi")
    plot!(line_scales, mssb, seriestype=:line, linewidth=linewidth, label="MSSB max lseg="*string(l_seg))
    plot!(line_scales, msmb, seriestype=:line, linestyle=:dashdot, linewidth=linewidth, label="MSMB max lseg="*string(l_seg))
    title!(title)

    return plot1
end


function plot_9bus_line_sweep(alg, dyn, mssb, msmb, line_scales, l_seg, size, linewidth, title, legendloc)

    plot1 = plot(line_scales, alg,xlabel="Line scale", ylabel="Max real λ != 0", label="Algebraic", legend=legendloc, seriestype=:line, linewidth=linewidth, size=size)
    plot!(line_scales, dyn, seriestype=:line, linewidth=linewidth, label="Dynpi")
    plot!(line_scales, mssb, seriestype=:line, linewidth=linewidth, label="MSSB max lseg="*string(l_seg))
    plot!(line_scales, msmb, seriestype=:line, linestyle=:dashdot, linewidth=linewidth, label="MSMB max lseg="*string(l_seg))
    title!(title)

    return plot1
end


function plot_2bus_eigenvalue_comparison(file_name::String, p1, p3, load_bus::String)

    sim_alg = ETL.build_2bus_sim_from_file(file_name, false, false, p1,load_bus)
    eigs_alg = small_signal_analysis(sim_alg).eigenvalues;
    sim_dyn = ETL.build_2bus_sim_from_file(file_name, true, false, p1, load_bus)
    eigs_dyn = small_signal_analysis(sim_dyn).eigenvalues;
    sim_mssb = ETL.build_2bus_sim_from_file(file_name, true, true, p1, load_bus)
    eigs_mssb = small_signal_analysis(sim_mssb).eigenvalues;
    sim_msmb = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus)
    eigs_msmb = small_signal_analysis(sim_msmb).eigenvalues;
    
    plot(real(eigs_alg), imag(eigs_alg), seriestype=:scatter, xlabel="Real", ylabel="Imag", label="Algebraic",legend = :outertopright)
    plot!(real(eigs_dyn), imag(eigs_dyn), seriestype=:scatter,label="Dynpi")
    plot!(real(eigs_mssb), imag(eigs_mssb), seriestype=:scatter,label="MSSB")
    plot!(real(eigs_msmb), imag(eigs_msmb), seriestype=:scatter,label="MSMB")
    title!("System eigs, p="*string(p1.p_load*p1.load_scale)*", q="*string(p1.q_load*p1.load_scale))

    
end
