cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Revise
using LaTeXStrings
using EffectsOfTLDynamics
using PlotlyJS
using Plots

using CSV
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

file_name = "../data/json_data/inv_v_machine.json"; # choose system 
line_dict = default_2_bus_line_dict

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 0.1,
)


# Choose parameters 
l = 100; # base line length 
l_seg = 10.0; # km 
line_scale = 1.0;
load_scale = 1.0;
load_bus = "BUS 2";

line_dict["BUS 1-BUS 2-i_1"] = l;
line_dict["BUS 1-BUS 2-i_1_static"] = l;

## Get line parameters 

M = 1;
factor_z = 1.0;
factor_y = 1.0;
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y);

M = 3;
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y);

# Calculate SIL 
V_nom = 230; # kV
Z_o = sqrt(z_km_ω_5_to_1_1/y_km_1);
SIL = V_nom ^2/Z_o;
p_load = real(SIL)/100;
q_load = imag(SIL)/100;

p3 = ETL.ExpParams(
    nothing, # N
    3, # M
    l, # gets overwritten
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
    line_scale, # line_scale 
    load_scale # load_scale
);

function run_line_sweep()
    # Line scale sweep 
    p1.load_scale = 1.0;
    p3.load_scale = 1.0;
    colorp = Plots.palette(:default); # color palette 

    p1.line_scale = 1.0;
    p3.line_scale = 1.0;
    sim = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus)
    eigvals = small_signal_analysis(sim).eigenvalues
    # Remove 0 eig 
    eigvals = eigvals[1:end-1]
    plot!(real(eigvals), imag(eigvals), seriestype=:scatter, color=colorp[4], label=L"$1.0$",ms=12, msw=1.0)
    p1.line_scale = 3.0;
    p3.line_scale = 3.0;
    sim = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus)
    eigvals = small_signal_analysis(sim).eigenvalues;
    # Remove 0 eig 
    eigvals = eigvals[1:end-1];
    plot!(real(eigvals), imag(eigvals), seriestype=:scatter, color=:thistle1, label=L"$3.0$",ms=10, msw=0.6)

    ylabel!(L"\mathrm{Im}\ [\lambda]")
    xlabel!(L"\mathrm{Re}\ [\lambda]")
    plot!(framestyle=:box)

    # ZOOM
    xlims!(-130,2)
    ylims!(-1e5,1e5)
    ylims!(-.7e5, .7e5)

    # FONT SIZING
    labelfontsize=25;
    ticksize = 16;
    legendsize = 16;
    plot!(legendtitle=L"Line\ scale")
    display(plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize, legendtitlesize=legendsize, legend=:bottomright))
end

# ARROWS
pt1 = (-.5,800)
pt2 = (.5,800)
Plots.plot(size=(1000,800))
plot!([0,0], [-7e4,7e4], linestyle=:dash, color=:black, primary=false)
plot!([pt1,pt2], arrow=true, color=:red, linewidth=2, primary=false)

run_line_sweep()
lens!([-5, 2], [-5000, 5000];
inset=(1, bbox(0.62, 0.03, 0.32, 0.32)), 
subplot=2, ticks=true, framestyle=:box,
lw=2, ls=:dot, lc=:orange)


figtitle = "../figures/Ruth/current/eig_plots/twobus_msmb_line_sweep_zoomed_inset"
savefig(figtitle*".svg")
savefig(figtitle*".png")


function run_load_sweep()

    # Load scale sweep
    p1.line_scale = 1.0;
    p3.line_scale = 1.0;

    p1.load_scale = 1.0;
    p3.load_scale = 1.0;
    sim = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus)
    eigvals = small_signal_analysis(sim).eigenvalues
    # Remove 0 eig 
    eigvals = eigvals[1:end-1]
    plot!(real(eigvals), imag(eigvals), seriestype=:scatter, color=colorp[4], label=L"$1.0$",ms=12, msw=1.0)
    p1.load_scale = 3.0;
    p3.load_scale = 3.0;
    sim = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus)
    eigvals = small_signal_analysis(sim).eigenvalues
    # Remove 0 eig 
    deleteat!(eigvals,length(eigvals)-2)
    plot!(real(eigvals), imag(eigvals), seriestype=:scatter, color=:thistle1, label=L"$3.0$",ms=10, msw=0.6)

    ylabel!(L"\mathrm{Im}\ [\lambda]")
    xlabel!(L"\mathrm{Re}\ [\lambda]")
    plot!(framestyle=:box)
    # ZOOM
    xlims!(-130,2)
    ylims!(-.62e5, .62e5)

    # FONT SIZING
    labelfontsize=25;
    ticksize = 16;
    legendsize = 16;
    plot!(legendtitle=L"Line\ scale")
    display(plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize, legendtitlesize=legendsize, legend=:bottomright))
end

# ARROWS
pt1 = (-.5,800)
pt2 = (.5,800)
Plots.plot(size=(1000,800))
plot!([0,0], [-7e4,7e4], linestyle=:dash, color=:black, primary=false)
plot!([pt1,pt2], arrow=true, color=:red, linewidth=2, primary=false)

run_load_sweep()
lens!([-5, 2], [-5000, 5000];
inset=(1, bbox(0.62, 0.03, 0.35, 0.35)), 
subplot=2, ticks=true, framestyle=:box,
lw=2, ls=:dot, lc=:orange)


figtitle = "../figures/Ruth/current/eig_plots/twobus_msmb_load_sweep_zoomed_inset"

savefig(figtitle*".svg")
savefig(figtitle*".png")

