



function pesgm_build_sim(file_name, alg_df, dyn_df, mssb_df, msmb_df, gfm_bus, gfl_bus, sm_bus)

    # Sample parameters 
    kpv = rand(kpv_nrel_range); # inner loop voltage gain 
    kiv = rand(kiv_nrel_range); # inner loop integral gain 
    kq = rand(gfm_kq_range); # outer loop reactive power control gain

    kpc = rand(kpc_range); # Inner current loop gain - GFM 
    kic = rand(kic_range); # Inner current loop gain - GFM

    load_scale = rand(load_scale_range)
    inv_share = rand(inverter_share_range)
    gfm_share = rand(gfm_share_range)
    
    # Algebraic 
    sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share, gfm_bus, gfl_bus);
    stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, false, false, p1, perturbation1,  gfm_bus, gfl_bus, sm_bus)
    push!(alg_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

    # Dynamic 
    sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share, gfm_bus, gfl_bus);
    stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, true, false, p1, perturbation1,  gfm_bus, gfl_bus, sm_bus)
    push!(dyn_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

    # # MSSB 
    # sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share);
    # stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, true, true, p1, perturbation1)
    # push!(mssb_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

    # # MSMB 
    # sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share);
    # stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, true, true, p3, perturbation1)
    # push!(msmb_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

end
