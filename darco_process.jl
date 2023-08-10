# Implements the process from D'Arco paper - finds bandwidth of line and converter eigs 
include("./all_sims.jl");
include("extra_functions.jl");
using LinearAlgebra
using DelimitedFiles

##### Define system and parameters #####

file_name = "test_sys.json"; # Inverter vs infinite bus 

t_max = 30.0;
dist = "CRC"; # control reference change (on inverter)

Z_c = 380; # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
abstol = 1e-13
reltol = 1e-10
maxiters = Int(1e10)

N = 1; # segs
l = 100; #km
p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters);

#### Create simulation and do small signal analysis on it ########

sim = get_ms_dyn_sim(file_name, t_max, dist, p)
sys = get_ms_dyn_sys(file_name, p);
ss = small_signal_analysis(sim)

#writedlm("test_psid_jacobian_export.csv", ss.reduced_jacobian, ',');
#writedlm("test_state_index_export.csv", ss.index, ',');

#summary_eigenvalues(ss)
pf = ss.participation_factors; # get all participation factors  as dict 
pf_df = summary_participation_factors(ss); 

# Get inverter 2 states
inv_name = "generator-102-1";
inv = get_component(DynamicInverter, sys, inv_name);
#inv_pf = pf[inv_name]; # get inverter participation factors 

inv_df = filter(r -> any(occursin.([inv_name], r.Name)), pf_df); # slice section of summary df that contains inverter states 

n_inv = zeros(1,size(pf_df,1)); # empty array for inverter interaction modes

for i = 1:size(pf_df,1) # for each mode
    p_i_subs = norm(inv_df[!,i+1], 1); # calculate norm of subsystem participation factors in that mode
    p_i = norm(pf_df[!,i+1], 1); #calculate norm of all participation factors in that mode
    n_inv[i] = p_i_subs/p_i; # divide 
end

# Get our line states... 
# In this case, these are just the remaining states bc the system is inverter vs infinite bus. 
# Might be harder for a larger system. Would need to iterate through all line names 

line_df = filter(r -> !any(occursin.([inv_name], r.Name)), pf_df);

n_line = zeros(1,size(pf_df,1)); # empty array for line interaction modes

for i = 1:size(pf_df,1) # for each mode
    p_i_subs = norm(line_df[!,i+1], 1); # calculate norm of subsystem participation factors in that mode
    p_i = norm(pf_df[!,i+1], 1); #calculate norm of all participation factors in that mode
    n_line[i] = p_i_subs/p_i; # divide 
end

# Return the modes that both the converter and line participate in

thr = 0.5; # threshold for participation

coupled = (n_inv.>χ) .&  (n_line.>χ)

# Find eig values for coupled eigs

coupled_eigs = ss.eigenvalues.*coupled'; # what is a smarter way to do this?
# Find states that participate most in these eigs
