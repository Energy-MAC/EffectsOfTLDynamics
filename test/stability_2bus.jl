cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
using EffectsOfTLDynamics

using CSV
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("eig_analysis_funcs.jl")

line_dict = default_2_bus_line_dict
foldername = "../figures/Ruth/current/twobus";
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

### Extract line data from files for M=1 and M=3
M = 1
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

M = 3
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

# Calculate SIL 
V_nom = 230 # kV
Z_o = sqrt(z_km_ω_5_to_1_1/y_km_1)
SIL = V_nom ^2/Z_o
p_load = real(SIL)/100
q_load = imag(SIL)/100

l_seg = 10.0; # km 

p1 = ETL.ExpParams(
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
    1.0, # line_scale - gets overwritten
    1.0 # load_scale - gets overwritten
);

p3 = ETL.ExpParams(
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
    1.0, # line_scale - gets overwritten
    1.0 # load_scale - gets overwritten
);

##------------------##

## INVERTER VS MACHINE 
file_name = "../data/json_data/inv_v_machine.json"; 

load_bus = "BUS 2"; # LOAD AT INVERTER BUS 
l_seg = 2.0;
# load_scale = 3.0;
# line_scales = collect(0.1:0.05:1.0);

load_scales = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
line_scales = [5.0:0.1:6.0, 3.0:0.1:4, 1.0:0.1:2.0, 0.5:0.1:1.0, 0.4:0.02:0.6, 0.3:0.05:0.5]

#alg, dyn, mssb, msmb = generate_2bus_line_sweep_data(file_name, p1, p3, load_scale, line_scales, l_seg, load_bus)
#plot_2bus_line_sweep(alg, dyn, mssb, msmb, line_scales, (800,600), 2, "TITLE", :top)


#load_scales = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]; # PARAMETER RANGES FOR BASE = SIL 
#line_scales = [8.0:0.1:9.0, 2.0:0.1:3.0, 0.6:0.1:1.2, 0.0:0.1:1.0, 0.0:0.1:1.0, 0.0:0.1:1.0];

alg, dyn, mssb, msmb = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, line_scales, l_seg, load_bus)

plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, (800,600), 2, (0.0,7.0), "Line scale (base=100km)", "Loading", "TITLE", true)


savepath = foldername*"/inv_mach_load_bus2_SIL_new.jld2"
@save savepath alg dyn mssb msmb load_scales line_scales p1 p3 load_bus l_seg 

## INVERTER VS MACHINE 
load_bus = "BUS 1"; # LOAD AT MACHINE BUS 
l_seg = 5.0;
load_scale = 1.5;
line_scales = collect(8:0.5:9.5);

alg, dyn, mssb, msmb = generate_2bus_line_sweep_data(file_name, p1, p3, load_scale, line_scales, l_seg, load_bus)
plot_2bus_line_sweep(alg, dyn, mssb, msmb, line_scales, (800,600), 2, "TITLE", :top)

l_seg = 5.0;
load_scales = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
line_scales = [6.5:0.1:7.5, 7.5:0.1:8.5, 8.5:0.1:9.5, 8.7:0.1:9.7, 8.0:0.1:9.0, 7.5:0.1:8.5]

alg, dyn, mssb, msmb = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, line_scales, l_seg, load_bus)

plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, (800,600), 2, (1.0,10.0), "Line scale (base=100km)", "Loading", "TITLE", true)

savepath = foldername*"/inv_mach_load_bus1_SIL_new.jld2"
@save savepath alg dyn mssb msmb load_scales line_scales p1 load_bus l_seg 


##------------------##

# INVERTER v INVERTER 

file_name = "../data/json_data/twobus_2inv.json"; # choose system 
load_bus = "BUS 2"; 

load_scales = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
line_scales = [9.0:0.1:10.0, 8.0:0.1:9.0, 5.0:0.1:6.0, 3.5:0.1:4.5, 2.5:0.1:3.5, 2.0:0.1:3.0];

alg, dyn, mssb, msmb = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, line_scales, l_seg, load_bus)

plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, (800,600), 2, (1.0,10.0), "Line scale (base=100km)", "Loading", "TITLE", true)

savepath = foldername*"/inv_inv_load_bus2_SIL_new.jld2"
@save savepath alg dyn mssb msmb load_scales line_scales p1 load_bus l_seg 

##------------------##

# MACHINE V MACHINE
file_name = "../data/json_data/mach_v_mach.json"; # choose system 
perturbation_params.crc_params = CRCParam(DynamicGenerator, "generator-102-1", :V_ref, 0.95) # have to update source type to be dynamic generator 
load_bus = "BUS 2"; 

l_seg = 5.0;
load_scales = [0.5, 1.0, 1.5, 2.0]
line_scales = [5.0:0.1:6.0, 4.0:0.1:5.0, 4.0:0.1:4.5, 0.1:0.1:0.3]

alg, dyn, mssb, msmb = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, line_scales, l_seg, load_bus)

plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, (800,600), 2, (0.0,7.0), "Line scale (base=100km)", "Loading", "TITLE", true)

savepath = foldername*"/mach_v_mach_load_bus2_SIL_new.jld2"
@save savepath alg dyn mssb msmb load_scales line_scales p1 load_bus l_seg 

# PLOT RESULTS 
cd(foldername)
ylims = (0.0,11.0);
xlims = (0.5,3.08);
ylabel = L"\mathrm{Line}\ \mathrm{scale}";
xlabel = L"\mathrm{Load}\ \mathrm{scale}";
figsize = (1000,800);
lw = 6; # line width 

# FONT SIZING
labelfontsize=40;
ticksize = 30;
legendsize = 24;

@load "inv_mach_load_bus2_SIL_new.jld2"
plot1 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, figsize, lw, ylims, xlims, ylabel, "", "", true)
plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize,legend=:topright)


figname = "inv_mach_l2"
savefig(figname*".svg")
savefig(figname*".png")

@load "inv_mach_load_bus1_SIL_new.jld2"
plot2 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales,figsize, lw, ylims, xlims, "", "", "", false)
plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize)


figname = "inv_mach_l1"
savefig(figname*".svg")
savefig(figname*".png")

@load "mach_v_mach_load_bus2_SIL_new.jld2"
plot3 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, figsize, lw, ylims, xlims, ylabel, xlabel, "", false)
plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize)


figname = "mach_mach_l2"
savefig(figname*".svg")
savefig(figname*".png")


@load "inv_inv_load_bus2_SIL_new.jld2"
plot4 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, figsize, lw, ylims, xlims, "", xlabel, "", false)
plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize)


figname = "inv_inv_l2"
savefig(figname*".svg")
savefig(figname*".png")


