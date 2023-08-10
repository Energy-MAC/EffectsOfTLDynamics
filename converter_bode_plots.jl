include("all_sims.jl")
include("extra_functions.jl")
using PowerSystems
using InfrastructureSystems
using InvertedIndices
using ControlSystems
using LinearAlgebra

sys = get_modified_twobus_sys(); # get inv vs inv system 

# Define sim parameters
t_max = 5.0;
tspan = (0.0, t_max);
dist = "CRC";
perturbation = choose_disturbance(sys, dist);

# Define line parameters 
Z_c = 380; # Ω
r_km = 0.05; # Ω/km
x_km = 0.488; # Ω/km
g_km = 0; # S/km
b_km = 3.371e-6; # S/km
abstol = 1e-13;
reltol = 1e-10;
maxiters = Int(1e10);

l = 50; # line length km

Nrange = [1,2,3,4,5,6,7,8,9,10];

#LsegRange = [10, 20, 25, 50, 100];
results = zeros(length(Nrange), 6);
j = 1;
for N = Nrange;

    ## - try to extract initial state values
    p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
    og_sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
    sys_ms = build_seg_model(og_sys_copy, p)
    sim = build_sim(sys_ms, tspan, perturbation, true);
    s = small_signal_analysis(sim);

    J = s.reduced_jacobian;

    state_index = collect(1:size(J)[1]); # Get list of all state indices 

    # iterate through 'devices' in model to get corresponding states for those devices. 
    d = get_state_index_map(s.index);

    # get indices for each susbsytem 

    line_states = get_subsystem_states!(["V_", "seg"], d)
    line_states_names = find_key_from_val!(d, line_states)

    inv1_states = get_subsystem_states!(["101"], d)
    inv1_states_names = find_key_from_val!(d,inv1_states)

    inv2_states = get_subsystem_states!(["102"], d)
    inv2_states_names = find_key_from_val!(d,inv2_states)

    # Get indices corresponding to the A and B matrices 
    inv1_A_i = inv1_states;
    inv1_B_i = setdiff(state_index, inv1_A_i);

    inv2_A_i = inv2_states;
    inv2_B_i = setdiff(state_index, inv2_A_i);

    line_A_i = line_states;
    line_B_i = setdiff(state_index, line_A_i);

    # Find which states need to be inputs for which subsystm - input_indices in global frame 
    line_input_indices, lineBframe = get_input_indices(J, line_A_i, line_B_i);
    inv1_input_indices, inv1Bframe = get_input_indices(J, inv1_A_i, inv1_B_i);
    inv2_input_indices, inv2Bframe = get_input_indices(J, inv2_A_i, inv2_B_i);

    # Find which states need to be outputs for each subsystem - output_indices in global frame 
    inv1_output_indices =  intersect!(vcat(inv2_input_indices, line_input_indices), inv1_states);
    inv2_output_indices = intersect!(vcat(line_input_indices, inv1_input_indices), inv2_states);
    line_output_indices = intersect!(vcat(inv1_input_indices, inv2_input_indices), line_states);

    # Construct subsystem models
    lineA = construct_A_matrix(J, line_A_i);
    lineB = construct_B_matrix(J, line_A_i, line_B_i, lineBframe);
    lineC = construct_C_matrix(line_A_i,line_output_indices);

    #line_inputs_names = find_key_from_val!(d, line_input_indices);
    #line_outputs_names = find_key_from_val!(d,line_output_indices);

    #Gline = ControlSystems.ss(lineA, lineB, lineC, 0);

    # Find line bandwidth - max imaginary eig part 
    print("Line subsys bandwidth: ")
    line_bw = maximum(abs.(imag(eigvals(lineA))))

    # # Construct inverter state space 
    inv1A = construct_A_matrix(J, inv1_A_i);
    inv2A = construct_A_matrix(J, inv2_A_i);

    # Find converter bandwidth 
    display("Inv1 subsys bandwidth: ")
    inv1_bw = maximum(abs.(imag(eigvals(inv1A))))

    display("Inv2 subsys bandwidth: ")
    inv2_bw = maximum(abs.(imag(eigvals(inv2A))))

    display("Combined system bandwidth")
    sys_bw = maximum(abs.(imag(s.eigenvalues)))

    # Do Darco proces of finding interaction modes. 

    pf = summary_participation_factors(s); # get all pf 
    # Need to find the indices of the pf_df that correspond to different subsystems.
    # These are different to indexing in state space.  

    line_df = get_subsys_df(["seg","V_"], pf)
    inv1_df = get_subsys_df(["101"], pf)
    inv2_df = get_subsys_df(["102"], pf)

    line_n = get_interaction_modes(pf, line_df)
    inv1_n = get_interaction_modes(pf, inv1_df)
    inv2_n = get_interaction_modes(pf, inv2_df)

    # Find modes that both inverters participate in
    thr = 0.05; # threshold for participation
    coupled = (inv1_n.>thr) .&  (inv2_n.>thr);

    # Find eig values for coupled eigs
    coupled_eigs = s.eigenvalues.*coupled'; # what is a smarter way to do this?
    count(real(coupled_eigs).<0) # see how many coupled eigs there are 
    print("Bandwidth of interacting inverter modes: ")
    interacting_bw = maximum(abs.(imag(coupled_eigs))) # find bandwidth of interacting modes 
    # Calculate the number of modes you can eliminate by using interacting_bw
    mode_savings = (length(s.eigenvalues) - count(real(coupled_eigs).<0))/length(s.eigenvalues);
    results[N,:] = [line_bw, inv1_bw, inv2_bw, sys_bw, interacting_bw, mode_savings]
    j += 1

end
results
xax = Nrange;

plot(xax, results[:,1], seriestype=:scatter, label="Line bandwidth")
plot!(xax, results[:,2], seriestype=:scatter, label="Inv bandwidth")
#plot!(xax, results[:,3], label="Inv2 bandwidth")
#plot!(xax, results[:,4], label="Sys bandwidth")
plot!(xax, results[:,5], seriestype=:scatter, label="Bandwidth from D'Arco")
xlabel!("N")
ylabel!("Bandwidth")
title!("Line length: "*string(l)*"km")

plot(xax, results[:,6]*100, label="% reduction in # states")
xlabel!("N")

function get_interaction_modes(pf, subsys_df)

    n = zeros(1,size(pf,1)); # empty array for line interaction modes
    for i = 1:size(pf,1) # for each mode
        p_i_subs = norm(subsys_df[!,i+1], 1); # calculate norm of subsystem participation factors in that mode
        p_i = norm(pf[!,i+1], 1); #calculate norm of all participation factors in that mode
        n[i] = p_i_subs/p_i; # divide 
    end
    return n
end

# inv1B = construct_B_matrix(J, inv1_A_i, inv1_B_i, inv1Bframe)
# inv1C = construct_C_matrix(inv1_A_i, inv1_output_indices)
# Ginv1 = ControlSystems.ss(inv1A, inv1B, inv1C, 0)


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

function find_key_from_val(x::Dict, target_val::Int);
    for (k,v) in x;
        if v==target_val;
            return k
            
        end
    end
end

function find_key_from_val!(x::Dict, target_val::Vector{Int});
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
