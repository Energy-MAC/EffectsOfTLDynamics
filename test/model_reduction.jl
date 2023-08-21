cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using EffectsOfTLDynamics

using InvertedIndices
using ControlSystems
using LinearAlgebra

# Define functions

function get_interaction_modes(pf, subsys_df)

    n = zeros(1,size(pf,1)); # empty array for line interaction modes
    for i = 1:size(pf,1) # for each mode
        p_i_subs = norm(subsys_df[!,i+1], 1); # calculate norm of subsystem participation factors in that mode
        p_i = norm(pf[!,i+1], 1); #calculate norm of all participation factors in that mode
        n[i] = p_i_subs/p_i; # divide 
    end
    return n
end

# Construct A, B, C matrices for each subsystem
function construct_C_matrix(Ai, indices)
    # C will be Ai * i dimension
    C = zeros(length(indices),size(Ai,1))
    for i in 1:size(indices,1);
        # Create a row with a 1 at the state index column 
        pos = findall(x->x==indices[i],Ai)
        C[i,pos[1]] = 1
    end
    return C
end

function construct_B_matrix(J,Ai, Bi, Bframe);
    B = J[Ai,Bi];
    B = B[:,Bframe]
    return B
end

function construct_A_matrix(J,Ai);
    A = J[Ai,Ai]
    return A
end

# Get B indices
function get_input_indices(J, Ai, Bi)
    B = J[Ai, Bi];
    input_states = findall(!iszero,B);
    input_indices_Bframe = []
    for i in input_states;
        push!(input_indices_Bframe, i[2])
    end
    input_indices = Bi[input_indices_Bframe]
    return unique(input_indices), unique(input_indices_Bframe)
end 

function get_state_index_map(x);
    d = Dict()
    for key in keys(x);
        for (state, i) in x[key];
            name = key*" "*String(state)
            d[name] = i;
        end
    end
    return d
end

function get_val!(x::Dict, target_val::Vector{Int});
    output = []
    for tval in target_val;
        for (k,v) in x;
            if v==tval;
                push!(output,k)
            end
        end
    end
    return output
end

function get_subsystem_states(x::String, d::Dict);
    # x is a string naming a subsystem 
    y = [];
    for (k, v) in d;
        if occursin(x, k);
            push!(y,v)
        end
    end
    return y 
end

function get_subsystem_states!(x::Vector{String}, d::Dict);
    # x is a vector of strings 
    y = [];
    for i in x;
        yi = get_subsystem_states(i, d)
        push!(y, yi)
    end
    return collect(Iterators.flatten(y)) # return as vector 

end

function get_subsys_df(x::Vector{String}, pf_df);
    y = []
    for i in x;
        for row in 1:size(pf_df, 1);
            state_name = pf_df[row,"Name"];
            if occursin(i,state_name)
                push!(y,row)
            end
        end
    end
    return pf_df[y,:]
end 




file_name = "../data/json_data/twobus_2inv.json"; # choose system 
sys = System(joinpath(pwd(), file_name));

# Define list of strings that can be used to uniquely identify each subsystem (e.g., 102, 101)
inv_ids = ["101", "102"]; 
line_ids = [["BUS 1-BUS 2", "V_"]]; # strings to identify line - need V_ to pick up voltage states. 

# Specify what subsystems are connected by specifying pairs of subsystem names - can be automated in future 
connections = [("101", "BUS 1-BUS 2"), ("102", "BUS 1-BUS 2")];

# Define sim parameters 
t_max = 2.0;
sim_p = SimParams(
    abstol = 1e13, # should this be -13? 
    reltol = 1e10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = t_max,
)

# Define line parameters 
Z_c = 380; # Ω
r_km = 0.05; # Ω/km
x_km = 0.488; # Ω/km
g_km = 0; # S/km
b_km = 3.371e-6; # S/km

## Commented out bc of issue with L shaped segments 

# l = 50; # line length km
# N = 3; # Number of segments. 
# M = 1; # Number of parallel branches 

# Perturbation type doesn't matter for this, just need to define something. 
#perturbation_type = "CRC"
#t_fault = 0.25
#default_crc_params = CRCParam(DynamicInverter, "generator-102-1", :V_ref, 0.95)
#perturbation_params = PerturbationParams(t_fault, nothing, nothing, default_crc_params, nothing, nothing, nothing)

#p = ExpParams(N, M, l, Z_c, r_km, x_km, g_km, b_km, sim_p, perturbation_type, perturbation_params)

og_sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
#sys_ms = build_seg_model!(og_sys_copy, p)
#dist = choose_disturbance(og_sys_copy, perturbation_type, p)
#sim = build_sim(sys_ms, (0.0,t_max), dist, true, p);

# Workaround while there's issues w L shape segments 
g = get_component(DynamicInjection, sys, "generator-102-1")
crc = ControlReferenceChange(1.0, g, :P_ref, 0.95)

# Make a pi line model 
sim = PSID.Simulation(
           MassMatrixModel, #Type of model used
           og_sys_copy, #system
           pwd(), #folder to output results
           (0.0, t_max), #time span
           crc, #Type of perturbation
           all_lines_dynamic = true
       )

ssig = small_signal_analysis(sim);

J = ssig.reduced_jacobian;

state_index = collect(1:size(J)[1]); # Get list of all state indices 

# iterate through 'devices' in model to get corresponding states for those devices. 
d = get_state_index_map(ssig.index);

# Get info on all the line subsystems 
for line_subsys in line_ids;
    line_states = get_subsystem_states!(line_subsys, d)
    display(get_val!(d, line_states))
    line_A_i = line_states;
    lineA = construct_A_matrix(J, line_A_i);
    line_bw = maximum(abs.(imag(eigvals(lineA))))
    display("Line bandwidth: "*string(line_bw))
end

# Get info on all the inverter subsystems 
for inv_subsys in inv_ids;
    inv_states = get_subsystem_states!([inv_subsys], d)
    display(get_val!(d, inv_states))
    inv_A_i = inv_states;
    invA = construct_A_matrix(J, inv_A_i);
    inv_bw = maximum(abs.(imag(eigvals(invA))))
    display("Inverter bandwidth: "*string(inv_bw))
end

display("Combined system bandwidth")
sys_bw = maximum(abs.(imag(ssig.eigenvalues)))

# Darco proces of finding interaction modes. 

pf = summary_participation_factors(ssig); # get all pf 
# Need to find the indices of the pf_df that correspond to different subsystems.
# These are different to indices in state space.  

subsys_ns = Dict(); # store vectors of overall participation by each subsys 
for subsys in union(inv_ids, line_ids);
    df = get_subsys_df([subsys], pf);
    n = get_interaction_modes(pf, df)
    display(n)
    subsys_ns[subsys] = n;
end

# Define threshold 
thr = 0.05; # threshold for participation

interacting_modes = [];
for pair in connections;
    s1 = subsys_ns[pair[1]];
    s2 = subsys_ns[pair[2]];
    coupled = (s1.>thr) .&  (s1.>thr);
    coupled_eigs = ssig.eigenvalues[vec(coupled)]; # what is a smarter way to do this?
    append!(interacting_modes, coupled_eigs)
end

interacting_modes = unique(interacting_modes);
print("Bandwidth of interacting inverter modes: ")
interacting_bw = maximum(abs.(imag(interacting_modes)))

