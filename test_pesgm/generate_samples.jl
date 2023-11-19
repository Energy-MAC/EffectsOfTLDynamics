using CSV
using DataFrames
using Random 
### BASELINE DATA ####
# Parameters that will give stable system 
col_names = ["line_scale", "load_scale", "Ta", "kd", "kq", "kpv", "kiv", "kpc", "kic", "inv_share", "gfm_share"]

function generate_default_darco_params_csv()
    kq = 0.2;
    kpv = 0.59 # d'arco 
    kiv = 739.0 # d'arco 
    kpc = 1.27
    kic = 14.3
    load_scale = 1.0;
    line_scale = 1.0;
    inv_share = 0.34;
    gfm_share = 0.77;
    Ta = 2.0;
    kd = 400.0;

    
    df = DataFrame([name => [] for name in col_names])
    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
    CSV.write("default_darco_params.csv", df)
end

function generate_nrel_param_csv()

    #### PARAMETER SET 1 #####
    gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)
    gfm_Ta_range = collect(2.0)
    gfm_kd_range = collect(400)

    # Inner loop voltage - for droop GFM 
    kpv_range = collect(0:0.003:0.01) # UNIFI
    kiv_range = collect(3:4:15) # UNIFI  

    # Inner loop current gains 
    kpc_range = collect(0.74:0.15:1.27) #? min=markovic max=d'arco
    kic_range = collect(1.19:4.0:14.3) #? min=markovic, max=d'arco

    line_scale_range = collect(1.0)
    load_scale_range = collect(0.5:0.3:1.5)
    inverter_share_range = collect(0.0:0.2:1.0)
    gfm_share_range = collect(0.0:0.2:1.0)

    df = DataFrame([name => [] for name in col_names])

    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    CSV.write("nrel_parameters.csv", df)

end

function generate_nrel_kpv_adjusted_params_csv()
    # Inner loop voltage P gains in D'Arco/Markovic range, not NREL 
    gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)
    gfm_Ta_range = collect(2.0)
    gfm_kd_range = collect(400)

    # Inner loop voltage -
    kpv_range = collect(0.5:0.03:0.6) # D'Arco range 
    kiv_range = collect(3:4:15) # UNIFI range 

    # Inner loop current gains 
    kpc_range = collect(0.74:0.15:1.27) #? min=markovic max=d'arco
    kic_range = collect(1.19:4.0:14.3) #? min=markovic, max=d'arco

    line_scale_range = collect(1.0)
    load_scale_range = collect(0.5:0.3:1.5)
    inverter_share_range = collect(0.0:0.2:1.0)
    gfm_share_range = collect(0.0:0.2:1.0)

    df = DataFrame([name => [] for name in col_names])

    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    CSV.write("nrel_parameters_kpv_adjusted.csv", df)
end

function generate_nrel_kpv_kiv_adjusted_params_csv()
    # Inner loop voltage P gains in D'Arco/Markovic range, not NREL 
    gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)
    gfm_Ta_range = collect(2.0)
    gfm_kd_range = collect(400)

    # Inner loop voltage 
    kpv_range = collect(0.5:0.03:0.6) # D'Arco range 
    kiv_range = collect(400:100:800) # D'Arco range 

    # Inner loop current gains 
    kpc_range = collect(0.74:0.15:1.27) #? min=markovic max=d'arco
    kic_range = collect(1.19:4.0:14.3) #? min=markovic, max=d'arco

    line_scale_range = collect(1.0)
    load_scale_range = collect(0.5:0.3:1.5)
    inverter_share_range = collect(0.0:0.2:1.0)
    gfm_share_range = collect(0.0:0.2:1.0)

    df = DataFrame([name => [] for name in col_names])

    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    CSV.write("nrel_parameters_kpv_kiv_adjusted.csv", df)
end

function generate_conditioned_nrel_param_csv()
    gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)
    gfm_Ta_range = collect(0.5:0.4:2.0) # (min=NREL, max=d'arco)
    gfm_kd_range = collect(100:100:400) # (min=NREL, max=d'arco)

    # Inner loop voltage loop 
    kpv_range = collect(0.5:0.03:0.6) # range centred on D'Arco value (0.59?)
    kiv_range = collect(400:100:800) # range centred on D'Arco value (736) 

    # Inner loop current gains 
    kpc_range = collect(0.74:0.15:1.27) # min=markovic max=d'arco
    kic_range = collect(1.19:4.0:14.3) # min=markovic, max=d'arco

    line_scale_range = collect(1.0)
    load_scale_range = collect(0.5:0.3:1.5)
    inverter_share_range = collect(0.0:0.2:1.0)
    gfm_share_range = collect(0.0:0.2:1.0)

    df = DataFrame([name => [] for name in col_names])

    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    CSV.write("conditioned_parameters_nrel_darco.csv", df)
end

function generate_case_param_csvs(inv_share, low_load_scale, high_load_scale, low_gfl, high_gfl, line_scale) 
    ### ------- #####
    gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)
    gfm_Ta_range = collect(0.5:0.4:2.0) # (min=NREL, max=d'arco)
    gfm_kd_range = collect(100:100:400) # (min=NREL, max=d'arco)

    # Inner loop voltage loop 
    kpv_range = collect(0.5:0.03:0.6) # range centred on D'Arco value (0.59?)
    kiv_range = collect(400:100:800) # range centred on D'Arco value (736) 

    # Inner loop current gains 
    kpc_range = collect(0.74:0.15:1.27) # min=markovic max=d'arco
    kic_range = collect(1.19:4.0:14.3) # min=markovic, max=d'arco

    line_scale_range = collect(line_scale)
    inverter_share_range = collect(inv_share)
    load_scale_range = collect(low_load_scale)
    gfm_share_range = collect(1.0-low_gfl); 

    df = DataFrame([name => [] for name in col_names])

    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    CSV.write("low_load_low_gfl.csv", df)

    # FIXED GFM/GFL/LOAD 
    load_scale_range = collect(low_load_scale)
    gfm_share_range = collect(1.0-high_gfl)

    df = DataFrame([name => [] for name in col_names])
    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    CSV.write("low_load_high_gfl.csv", df)

    load_scale_range = collect(high_load_scale)
    gfm_share_range = collect(1.0-low_gfl)

    df = DataFrame([name => [] for name in col_names])
    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    CSV.write("nom_load_low_gfl.csv", df)

    load_scale_range = collect(high_load_scale)
    gfm_share_range = collect(1.0-high_gfl)

    df = DataFrame([name => [] for name in col_names])
    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    CSV.write("nom_load_high_gfl.csv", df)
end

function generate_single_case_params(inv_share, fixed_load_scale, gfl_share, line_scale) 
    ### ------- #####
    gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)
    gfm_Ta_range = collect(0.5:0.4:2.0) # (min=NREL, max=d'arco)
    gfm_kd_range = collect(100:100:400) # (min=NREL, max=d'arco)

    # Inner loop voltage loop 
    kpv_range = collect(0.5:0.03:0.6) # range centred on D'Arco value (0.59?)
    kiv_range = collect(400:100:800) # range centred on D'Arco value (736) 

    # Inner loop current gains 
    kpc_range = collect(0.74:0.15:1.27) # min=markovic max=d'arco
    kic_range = collect(1.19:4.0:14.3) # min=markovic, max=d'arco

    line_scale_range = collect(line_scale)
    inverter_share_range = collect(inv_share)
    load_scale_range = collect(fixed_load_scale)
    gfm_share_range = collect(1.0-gfl_share); 

    df = DataFrame([name => [] for name in col_names])

    for line_scale = line_scale_range;
        for kq = gfm_kq_range;
            for kpv = kpv_range;
                for kiv = kiv_range;
                    for kpc = kpc_range;
                        for kic = kic_range;
                            for load_scale = load_scale_range;
                                for inv_share = inverter_share_range;
                                    for gfm_share = gfm_share_range;
                                        for Ta = gfm_Ta_range;
                                            for kd = gfm_kd_range;
                                                # Check conditions 
                                                c1 = kic/kpc >= 10 # integral gain 
                                                c2 = kiv/kpv >= 10
                                                c3 = kpc > kpv 
                                                #c4 = kpv > kq #? unsure about this one 

                                                if c1 & c2 & c3
                                                    push!(df, [line_scale, load_scale, Ta, kd, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return df 
end

function generate_single_case_gain_params() 
    ### ------- #####
    gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)
    gfm_Ta_range = collect(0.5:0.4:2.0) # (min=NREL, max=d'arco)
    gfm_kd_range = collect(100:100:400) # (min=NREL, max=d'arco)

    # Inner loop voltage loop 
    kpv_range = collect(0.5:0.03:0.6) # range centred on D'Arco value (0.59?)
    kiv_range = collect(400:100:800) # range centred on D'Arco value (736) 

    # Inner loop current gains 
    kpc_range = collect(0.74:0.15:1.27) # min=markovic max=d'arco
    kic_range = collect(1.19:4.0:14.3) # min=markovic, max=d'arco

    col_names = ["Ta", "kd", "kq", "kpv", "kiv", "kpc", "kic"]

    df = DataFrame([name => [] for name in col_names])

    for kq = gfm_kq_range;
        for kpv = kpv_range;
            for kiv = kiv_range;
                for kpc = kpc_range;
                    for kic = kic_range;
                        for Ta = gfm_Ta_range;
                            for kd = gfm_kd_range;
                                # Check conditions 
                                c1 = kic/kpc >= 10 # integral gain 
                                c2 = kiv/kpv >= 10
                                c3 = kpc > kpv 
                                #c4 = kpv > kq #? unsure about this one 

                                if c1 & c2 & c3
                                    push!(df, [Ta, kd, kq, kpv, kiv, kpc, kic])
                                end 
                            end
                        end
                    end
                end
            end
        end
    end

    return df 
end

function get_random_sample(df, N)
    samples_shuffled = shuffle(MersenneTwister(123), df);
    samples_cut = samples_shuffled[1:N,:]
    return samples_cut
end 


