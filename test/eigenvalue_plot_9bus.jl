cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Revise
using LaTeXStrings
using Plots
using EffectsOfTLDynamics
#using PlotlyJS

using CSV
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

file_name = "../data/json_data/9bus_slackless.json"; # choose system 
line_dict = default_9_bus_line_dict

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

t_fault = 0.25
perturbation_type = "BranchTrip"
perturbation_params = get_default_perturbation(t_fault, perturbation_type)
perturbation_params.branch_trip_params = BTParam("Bus 7-Bus 5-i_1")

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

## Get line parameters 
M = 1
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

M = 3
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

p1 = ETL.ExpParams(
    nothing, # N 
    1, # M
    l, # gets overwritten
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
    0.0,
    0.0,
    line_scale, # line_scale 
    load_scale # load_scale 
);

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
    0.0,
    0.0,
    line_scale, # line_scale 
    load_scale # load_scale
);


sim_alg = ETL.build_9bus_sim_from_file(file_name, false, false, p1)
eigs_alg = small_signal_analysis(sim_alg).eigenvalues;
sim_dyn = ETL.build_9bus_sim_from_file(file_name, true, false, p1)
eigs_dyn = small_signal_analysis(sim_dyn).eigenvalues;
sim_mssb = ETL.build_9bus_sim_from_file(file_name, true, true, p1)
eigs_mssb = small_signal_analysis(sim_mssb).eigenvalues;
sim_msmb = ETL.build_9bus_sim_from_file(file_name, true, true, p3)
eigs_msmb = small_signal_analysis(sim_msmb).eigenvalues;



colorp = palette(:default); # color palette 
st = :scatter # plot style 
ms = 7; # marker size 
msw = 0.5; # marker size width 


plot1 = Plots.plot(real(eigs_alg), imag(eigs_alg), seriestype=st, ylabel="Imag", label=L"\mathrm{statpi}",legend = :outertopright, color=colorp[1], ms=ms, msw=msw, size=(1000,800))
plot!(real(eigs_dyn), imag(eigs_dyn), st=st,label=L"\mathrm{dynpi}", color=colorp[2], ms=ms,msw=msw)
plot!(real(eigs_mssb), imag(eigs_mssb), st=st,label=L"\mathrm{MSSB}", color=colorp[3], ms=ms, msw=msw)
plot!(real(eigs_msmb), imag(eigs_msmb), st=st,label=L"\mathrm{MSMB}", color=colorp[4], ms=ms, msw=msw)

#plot!(real(eigs_alg), imag(eigs_alg), seriestype=st, ylabel="Imag", label=L"\mathrm{statpi}",color=colorp[1])
#plot!(real(eigs_dyn), imag(eigs_dyn), st=st,label=L"\mathrm{dynpi}", color=colorp[2], ms=ms,msw=msw)

xlims!(-300,10)
ylims!(-1e5,1e5)
ylabel!(L"Imag")
xlabel!(L"Real")
title!(L"\mathrm{Load\ scale: }"*L"$%$load_scale,\ \mathrm{Line\ scale: }"*L"%$line_scale$")
#title!(L"Load scale: "*string(load_scale)*", line scale: "*string(line_scale))


#savefig("../figures/Ruth/current/ninebus_eigs_extra_zoom_2.png")


pf_df = summary_participation_factors(small_signal_analysis(sim_dyn))

pf_df[1:22, end-9:end]#.>0.01

# Filter by mode/eig 
filter(row -> row.λ_34 > 0.1, pf_df)[!,[:Name, :λ_34]]

# Filter by line states and look at their participation in least stable 5 eigs 
filter(row -> row.Name == "V_2 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_2 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "V_1 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_1 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_R", pf_df)[!,end-9:end]
filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_I", pf_df)[!,end-9:end]


pf_df = summary_participation_factors(small_signal_analysis(sim_msmb))

# Filter by line states and look at their participation in least stable 5 eigs 
# Line states are first 24 (?) states 
pf_df[1:82, end-9:end].>0.01


filter(row -> row.Name == "V_2 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_2 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "V_1 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_1 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_R", pf_df)[!,end-9:end]
filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_I", pf_df)[!,end-9:end]

