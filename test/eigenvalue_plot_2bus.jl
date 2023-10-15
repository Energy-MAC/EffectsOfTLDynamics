cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Revise
using LaTeXStrings
using EffectsOfTLDynamics
#using PlotlyJS
using Plots
using JLD2

using CSV
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

function draw_oval_with_edge(x, y, a, b, edgecolor, edgestrokealpha, edgestrokewidth)
    n_points = 100  # Number of points to approximate the oval shape
    θ = LinRange(0, 2π, n_points)  # Create an array of angles

    # Calculate the x and y coordinates of points on the oval
    x_points = a * cos.(θ) .+ x
    y_points = b * sin.(θ) .+ y

    # Plot the oval's edge
    plot!(x_points, y_points, linecolor=edgecolor, linealpha=edgestrokealpha, linewidth=edgestrokewidth, line=:path, primary=false)
end

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
    p_load,
    q_load,
    line_scale, # line_scale 
    load_scale # load_scale
);


sim_alg = ETL.build_2bus_sim_from_file(file_name, false, false, p1,load_bus)
eigs_alg = small_signal_analysis(sim_alg).eigenvalues;
sim_dyn = ETL.build_2bus_sim_from_file(file_name, true, false, p1, load_bus)
eigs_dyn = small_signal_analysis(sim_dyn).eigenvalues;
sim_mssb = ETL.build_2bus_sim_from_file(file_name, true, true, p1, load_bus)
eigs_mssb = small_signal_analysis(sim_mssb).eigenvalues;
sim_msmb = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus)
eigs_msmb = small_signal_analysis(sim_msmb).eigenvalues;

# # Save data 
# foldername = "../figures/Ruth/current/eig_plots"
# savepath = foldername*"/load1_line1_data.jld2"
# @save savepath sim_alg eigs_alg sim_dyn eigs_dyn sim_mssb eigs_mssb sim_msmb eigs_msmb p1 p3 load_bus 

# cd(foldername)
# @load "load1_line1_data.jld2"

colorp = Plots.palette(:default); # color palette 
st = :scatter # plot style 
ms = 10; # marker size 
msw = 1.0; # marker size width 

plot1 = Plots.plot(real(eigs_msmb), imag(eigs_msmb), st=st,label=L"\mathrm{MSMB}", color=colorp[4], ms=ms, msw=msw)
plot!(real(eigs_mssb), imag(eigs_mssb), st=st,label=L"\mathrm{MSSB}", color=colorp[3], ms=ms, msw=msw)
plot!(real(eigs_dyn), imag(eigs_dyn), st=st,label=L"\mathrm{dynpi}", color=colorp[2], ms=ms,msw=msw)
plot!(real(eigs_alg), imag(eigs_alg), seriestype=st, label=L"\mathrm{statpi}", color=colorp[1], ms=ms, msw=msw, size=(1000,800))

ylabel!(L"\mathrm{Im}\ [\lambda]")
xlabel!(L"\mathrm{Re}\ [\lambda]")

# ZOOM
xlims!(-200,10)
ylims!(-.72e5, .72e5)

# FONT SIZING
labelfontsize=25;
ticksize = 16;
legendsize = 16;
plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize,legend=:left)
plot!(framestyle=:box)

# MAKE ANNOTATIONS
draw_oval_with_edge(-12, 0, 20.0, 1e4, :red, 0.7, 2) 
annotate!(-15, -1e4, (L"GFM/SM\ states", 16, :red, :top))

draw_oval_with_edge(-95, 3.8e4, 100.0, 3.3e4, :red, 0.7, 2) 
annotate!(-85, 6.8e4, (L"Line\ states", 16, :red, :top))


#title!(L"\mathrm{Load\ scale: }"*L"$%$load_scale,\ \mathrm{Line\ scale: }"*L"%$line_scale$")
#title!(L"Load scale: "*string(load_scale)*", line scale: "*string(line_scale))
cd(@__DIR__)
figname = "../figures/Ruth/current/eig_plots/twobus_all_models_zoomed_3";
savefig(figname*".svg")
savefig(figname*".png")


# Look at dominant participation in eigs, for MSMB
eig_sum = summary_eigenvalues(small_signal_analysis(sim_msmb))
pf_df = summary_participation_factors(small_signal_analysis(sim_msmb))
pf_df[findall(pf_df."λ_1" .> 1e-5),["Name","λ_1"]]
pf_df[findall(pf_df."λ_2" .> 1e-5),["Name","λ_2"]]
pf_df[findall(pf_df."λ_3" .> 1e-5),["Name","λ_3"]]
pf_df[findall(pf_df."λ_4" .> 1e-5),["Name","λ_4"]]
pf_df[findall(pf_df."λ_5" .> 1e-5),["Name","λ_5"]]
pf_df[findall(pf_df."λ_6" .> 1e-5),["Name","λ_6"]]
pf_df[findall(pf_df."λ_7" .> 1e-5),["Name","λ_7"]]
pf_df[findall(pf_df."λ_8" .> 1e-5),["Name","λ_8"]]


pf_df[findall(pf_df."λ_113" .> 1e-5),["Name","λ_113"]]
pf_df[findall(pf_df."λ_112" .> 1e-5),["Name","λ_112"]]
pf_df[findall(pf_df."λ_8" .> 1e-5),["Name","λ_8"]]




# Sort by frequency 
sort!(eig_sum,["Freq [Hz]"], rev=true)[1:90,:]

# Sort by imag part > 1e3 
sort!(eig_sum, ["Imag. Part"])[1:50,:]






# PARTICIPATION FACTORS 
pf_df = summary_participation_factors(small_signal_analysis(sim_dyn))

# Filter by mode/eig 
filter(row -> row.λ_34 > 0.1, pf_df)[!,[:Name, :λ_34]]

# Filter by line states and look at their participation in least stable 5 eigs 
filter(row -> row.Name == "V_2 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_2 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "V_1 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_1 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_R", pf_df)[!,end-9:end]
filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_I", pf_df)[!,end-9:end]




# Filter by line states and look at their participation in least stable 5 eigs 
# Line states are first 24 (?) states 
pf_df[1:82, end-9:end].>0.01


filter(row -> row.Name == "V_2 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_2 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "V_1 R", pf_df)[!,end-9:end]
filter(row -> row.Name == "V_1 I", pf_df)[!,end-9:end]

filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_R", pf_df)[!,end-9:end]
filter(row -> row.Name == "BUS 1-BUS 2-i_1 Il_I", pf_df)[!,end-9:end]

