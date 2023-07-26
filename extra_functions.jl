# Define function to build simulation using multi-seg dynamic lines
function get_ms_dyn_sim(file_name, t_max, dist, p::ExpParams)
    sys = System(joinpath(pwd(), file_name));
    tspan = (0.0, t_max)
    perturbation = choose_disturbance(sys, dist)
    sys = build_seg_model(sys, p)
    sim = build_sim(sys, tspan, perturbation, true);
#    show_states_initial_value(sim);
    return sim
end


function get_ms_dyn_sys(file_name, p::ExpParams)
    sys = System(joinpath(pwd(), file_name));
    sys = build_seg_model(sys, p)
    return sys
end
