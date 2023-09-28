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
Z_c = sqrt(z_km_1[1]/y_km_1)
Vnom = 230; # kV 
SIL = Vnom^2/Z_c;
SIL_P = real(SIL)
SIL_Q = imag(SIL)

p_load = SIL_P/100 # convert to p.u. 
q_load = SIL_Q/100

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
# LOAD AT INVERTER BUS 
load_bus = "BUS 2";


# # DEFINE RANGES -- FOR BASE P=1.0, Q=0.25
# load_scales = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
# lRange = [1030:5:1060, 920:5:970, 700:5:760, 500:5:550, 360:5:380, 230:5:280, 160:5:220];
# # GENERATE DATA  
# alg, dyn, mssb, msmb, h = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, lRange, l_seg, load_bus)
# # SAVE DATA 
# foldername = "../figures/Ruth/twobus"
# savepath = foldername*"/inv_mach_load_bus2.jld2"
# @save savepath alg dyn mssb msmb h load_scales lRange p1 p3 load_bus l_seg 

# PARAMETER RANGES FOR BASE = SIL 
load_scales = [0.5, 0.75, 1.0, 1.25, 1.5]
lRange = [800:10:900,450:10:550, 200:10:300, 125:10:175, 60:10:120]
alg, dyn, mssb, msmb, h = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, lRange, l_seg, load_bus)

savepath = foldername*"/inv_mach_load_bus2_SIL.jld2"
@save savepath alg dyn mssb msmb h load_scales lRange p1 p3 load_bus l_seg 

## INVERTER VS MACHINE 

# LOAD AT MACHINE BUS 
load_bus = "BUS 1";
# # FOR BASE P=1.0, Q=0.25
# load_scales = [2.0, 2.5, 3.0, 3.5, 4.0];
# lRange = [720:5:800, 575:2:585, 440:5:460, 360:5:375, 290:5:310];
# alg, dyn, mssb, msmb, h = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, lRange, l_seg, load_bus)
# savepath = foldername*"/inv_mach_load_bus1.jld2"
# @save savepath2 alg dyn mssb msmb h load_scales lRange p1 p3 load_bus l_seg 

# PARAMETER RANGES FOR BASE = SIL 
load_scales = [0.5, 1.0, 1.5, 2.0];
lRange = [1175:10:1275, 1200:10:1300, 1100:10:1200, 800:10:900];
alg, dyn, mssb, msmb, h = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, lRange, l_seg, load_bus)
# h
# plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, p1,(800,600), 3)
savepath = foldername*"/inv_mach_load_bus1_SIL.jld2"
@save savepath alg dyn mssb msmb h load_scales lRange p1 load_bus l_seg 


##------------------##

# INVERTER v INVERTER 
file_name = "../data/json_data/twobus_2inv.json"; # choose system 
load_bus = "BUS 1"; 
load_scale = 1.0;
lRange = 480:10:530;
alg, dyn, mssb, msmb, h = generate_2bus_line_sweep_data(file_name, p1, p3, load_scale, lRange, l_seg, load_bus)
h


## INV V INV DEBUG 

l_seg = 5.0;
line_scales = Float64.(collect(1.0:1.0:8.0));
load_scale = 1.0;
alg1, dyn1, mssb1, msmb1 = generate_2bus_line_sweep_data(file_name, p1, p3, load_scale, line_scales, l_seg, "BUS 1");
plot1 = plot_2bus_line_sweep(alg1, dyn1, mssb1, msmb1, line_scales, l_seg, (800,600), 2, "LOAD @ BUS 1", :top);

alg2, dyn2, mssb2, msmb2 = generate_2bus_line_sweep_data(file_name, p1, p3, load_scale, line_scales, l_seg, "BUS 2");
plot2 = plot_2bus_line_sweep(alg2, dyn2, mssb2, msmb2, line_scales, l_seg, (800,600), 2, "LOAD @ BUS 2", :right)

l = @layout [a;b];
plot(plot1, plot2, layout=l)









# load_scales = [0.5, 1.0, 1.5, 2.0]
# lRange = [.., 480:5:530]
# alg, dyn, mssb, msmb, h = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, lRange, l_seg, load_bus)
# savepath = foldername*"/inv_inv_SIL.jld2"
# @save savepath alg dyn mssb msmb h load_scales lRange p1 load_bus l_seg 

##------------------##

# MACHINE V MACHINE
file_name = "../data/json_data/mach_v_mach.json"; # choose system 
perturbation_params.crc_params = CRCParam(DynamicGenerator, "generator-102-1", :V_ref, 0.95) # have to update source type to be dynamic generator 
load_bus = "BUS 1"; 
# FOR NON SIL PARAMS 
# load_scales = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
# lRange = [1030:5:1060, 920:5:970, 700:5:760, 500:5:550, 360:5:380, 230:5:280, 160:5:220];
# perturbation_params.crc_params = CRCParam(DynamicGenerator, "generator-102-1", :V_ref, 0.95) 
# alg, dyn, mssb, msmb, h = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, lRange, l_seg, load_bus)
# savepath = foldername*"/mach_mach.jld2"
# @save savepath alg dyn mssb msmb h load_scales lRange p1 load_bus l_seg 

# FOR SIL PARAMS 
# load_scale = 2.0;
# lRange = 200:50:400;
# #perturbation_params.crc_params = CRCParam(DynamicGenerator, "generator-102-1", :V_ref, 0.95)
# alg, dyn, mssb, msmb, h = generate_2bus_line_sweep_data(file_name, p1, p3, load_scale, lRange, l_seg, load_bus);
# h

load_scales = [0.5, 1.0, 1.5, 2.0]
lRange = [980:10:1050, 800:10:850, 620:10:670, 290:10:340]
alg, dyn, mssb, msmb = generate_2bus_loading_v_boundary_data(file_name, p1, p3, load_scales, lRange, l_seg, load_bus)

plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, p1,(800,600), 3)

savepath = foldername*"/mach_mach_SIL.jld2";
@save savepath alg dyn mssb msmb h load_scales lRange p1 load_bus l_seg 


## MACH v MACH DEBUG 

l_seg = 5.0;
line_scales = Float64.(collect(1.0:1.0:8.0));
alg1, dyn1, mssb1, msmb1 = generate_2bus_line_sweep_data(file_name, p1, p3, 1.0, line_scales, l_seg, "BUS 1");
plot1 = plot_2bus_line_sweep(alg1, dyn1, mssb1, msmb1, line_scales, l_seg, (800,600), 2, "LOAD @ BUS 1", :top);

alg2, dyn2, mssb2, msmb2 = generate_2bus_line_sweep_data(file_name, p1, p3, 1.0, line_scales, l_seg, "BUS 2");
plot2 = plot_2bus_line_sweep(alg2, dyn2, mssb2, msmb2, line_scales, l_seg, (800,600), 2, "LOAD @ BUS 2", :left);

l = @layout [a;b];
plot(plot1, plot2, layout=l)





# Produce 2x2 plot 

l = @layout [a b; c d]

cd(foldername)
ylims = (50,1500);
@load "inv_mach_load_bus2_SIL.jld2"
plot1 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, p1,(800,600), 3, ylims, "Line length (km)", "", "Inv v mach, l2", true)

@load "inv_mach_load_bus1_SIL.jld2"
plot2 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, p1,(800,600), 3, ylims, "", "", "Inv v mach, l1", false)

@load "mach_mach_SIL.jld2"

plot3 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, p1,(800,600), 3, ylims, "Line length (km)", "Load scale","Mach v mach", false)
plot4 = plot_loading_v_boundary(alg, dyn, mssb, msmb, load_scales, l_seg, p1,(800,600), 3, ylims, "", "Load scale", "Mach v mach", false)

plot(plot1, plot2, plot3, plot4, layout = l)
