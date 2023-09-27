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
l = 200;
l_seg = 10.0; # km 
line_scale = 1.0;
load_scale = 1.0;
load_bus = "BUS 2"
# Update dict 
line_dict["BUS 1-BUS 2-i_1"] = l;
line_dict["BUS 1-BUS 2-i_1_static"] = l;

## Get line parameters 

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

p_load = SIL_P/100;
q_load = SIL_Q/100;

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
    p_load,
    q_load,
    1.0, # line_scale 
    1.0 # load_scale 
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
    p_load,
    q_load,
    1.0, # line_scale 
    1.0 # load_scale
);


sim_alg = ETL.build_2bus_sim_from_file(file_name, false, false, p1,load_bus)
eigs_alg = small_signal_analysis(sim_alg).eigenvalues;
sim_dyn = ETL.build_2bus_sim_from_file(file_name, true, false, p1, load_bus)
eigs_dyn = small_signal_analysis(sim_dyn).eigenvalues;
sim_mssb = ETL.build_2bus_sim_from_file(file_name, true, true, p1, load_bus)
eigs_mssb = small_signal_analysis(sim_mssb).eigenvalues;
sim_msmb = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus)
eigs_msmb = small_signal_analysis(sim_msmb).eigenvalues;

colorp = palette(:default);
st = :scatter
ms = 8;
msw = 0;

plot1 = plot(real(eigs_alg), imag(eigs_alg), seriestype=seriestype, ylabel="Imag", label="Algebraic",legend = :outertopright, color=colorp[1], ms=ms, msw=msw, size=(1000,800));
plot!(real(eigs_dyn), imag(eigs_dyn), st=st,label="Dynpi", color=colorp[2], ms=ms-2,msw=msw)
plot!(real(eigs_mssb), imag(eigs_mssb), st=st,label="MSSB", color=colorp[3], ms=ms-4, msw=msw)
plot!(real(eigs_msmb), imag(eigs_msmb), st=st,label="MSMB", color=colorp[4], ms=ms-6, msw=msw)


plot2 = plot(real(eigs_mssb), imag(eigs_mssb), st=st,label="MSSB", color=colorp[3], legend = :outertopright, ms=ms,msw=msw);
plot!(real(eigs_msmb), imag(eigs_msmb), st=st,label="MSMB", color=colorp[4], ms=ms-2,msw=msw);

plot3 = plot(real(eigs_alg), imag(eigs_alg), st=st, xlabel="Real", ylabel="Imag", label="Algebraic",legend = :outertopright, color=colorp[1], ms=ms,msw=msw);
plot!(real(eigs_dyn), imag(eigs_dyn), st=st, label="Dynpi", color=colorp[2], ms=ms-2,msw=msw);

l = @layout [a;b;c];
plot(plot1, plot2, plot3, layout=l, size=(1000,800))
xlims!(-100,0.5)



# #ylims!(-4,4)

# title!("System eigs, p="*string(round(p_load*load_scale, digits=2))*", q="*string(round(q_load*load_scale, digits=2)))



# # Plot eigs with participation from different subsystems (e.g., line)
# # Find line states

#summary_eigenvalues(small_signal_analysis(sim_alg))

pf_df = summary_participation_factors(small_signal_analysis(sim_dyn))

# Filter by mode/eig 
filter(row -> row.λ_36 > 0.1, pf_df)[!,[:Name, :λ_36]]

# Filter by line states and look at their participation in least stable 5 eigs 
filter(row -> row.Name == "V_2 R", pf_df)[!,end-4:end]
filter(row -> row.Name == "V_2 I", pf_df)[!,end-4:end]

filter(row -> row.Name == "V_1 R", pf_df)[!,end-4:end]
filter(row -> row.Name == "V_1 I", pf_df)[!,end-4:end]

filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_R", pf_df)[!,end-4:end]
filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_I", pf_df)[!,end-4:end]



pf_df[!,pf_df[!,:λ_1]==0.0]
summary_participation_factors(small_signal_analysis(sim_alg))[:,[:Name, :λ_30]]