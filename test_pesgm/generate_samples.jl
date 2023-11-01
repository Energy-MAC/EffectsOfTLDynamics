
using CSV
using DataFrames

# Generate different sets of parameter samples 

### BASELINE DATA ####
# Parameters that will give stable system 
kq = 0.2;
kpv = 0.59 # d'arco 
kiv = 739.0 # d'arco 
kpc = 1.27
kic = 14.3
load_scale = 1.0;
inv_share = 0.34;
gfm_share = 0.77;

col_names = ["load_scale", "kq", "kpv", "kiv", "kpc", "kic", "inv_share", "gfm_share"]
df = DataFrame([name => [] for name in col_names])

push!(df, [load_scale, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])

CSV.write("default_darco_params.csv", df)

# ### BASELINE DATA ####
# # Parameters that will give stable system 
# kq = 0.2;
# kpv = 0.59 # d'arco 
# kiv = 739.0 # d'arco 
# kpc = 1.27
# kic = 14.3
# load_scale = 1.0;
# inv_share = 0.8;
# gfm_share = 0.77;

# col_names = ["load_scale", "kq", "kpv", "kiv", "kpc", "kic", "inv_share", "gfm_share"]
# df = DataFrame([name => [] for name in col_names])

# push!(df, [load_scale, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])

# CSV.write("test_params1.csv", df)



#### PARAMETER SET 1 #####
gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)

# Inner loop voltage - for droop GFM 
kpv_range = collect(0:0.003:0.01) # UNIFI
kiv_range = collect(3:4:15) # UNIFI  

# Inner loop current gains 
kpc_range = collect(0.74:0.15:1.27) #? min=markovic max=d'arco
kic_range = collect(1.19:4.0:14.3) #? min=markovic, max=d'arco

load_scale_range = collect(1.0)
inverter_share_range = collect(0.0:0.2:1.0)
gfm_share_range = collect(0.0:0.2:1.0)

col_names = ["load_scale", "kq", "kpv", "kiv", "kpc", "kic", "inv_share", "gfm_share"]
df = DataFrame([name => [] for name in col_names])

for kq = gfm_kq_range;
    for kpv = kpv_range;
        for kiv = kiv_range;
            for kpc = kpc_range;
                for kic = kic_range;
                    for load_scale = load_scale_range;
                        for inv_share = inverter_share_range;
                            for gfm_share = gfm_share_range;
                                push!(df, [load_scale, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                            end
                        end
                    end
                end
            end
        end
    end
end

CSV.write("nrel_parameters.csv", df)

### PARAMETER SET 2 #### 
# Inner loop voltage P gains in D'Arco/Markovic range, not NREL 
gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)

# Inner loop voltage -
kpv_range = collect(0.5:0.03:0.6) # D'Arco range 
kiv_range = collect(3:4:15) # UNIFI range 

# Inner loop current gains 
kpc_range = collect(0.74:0.15:1.27) #? min=markovic max=d'arco
kic_range = collect(1.19:4.0:14.3) #? min=markovic, max=d'arco

load_scale_range = collect(1.0)
inverter_share_range = collect(0.0:0.2:1.0)
gfm_share_range = collect(0.0:0.2:1.0)

col_names = ["load_scale", "kq", "kpv", "kiv", "kpc", "kic", "inv_share", "gfm_share"]
df = DataFrame([name => [] for name in col_names])

for kq = gfm_kq_range;
    for kpv = kpv_range;
        for kiv = kiv_range;
            for kpc = kpc_range;
                for kic = kic_range;
                    for load_scale = load_scale_range;
                        for inv_share = inverter_share_range;
                            for gfm_share = gfm_share_range;
                                push!(df, [load_scale, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                            end
                        end
                    end
                end
            end
        end
    end
end

CSV.write("nrel_parameters_kpv_adjusted.csv", df)


### ------- #####
# Inner loop voltage P gains in D'Arco/Markovic range, not NREL 
gfm_kq_range = collect(0.05:0.05:0.2) # (min=NREL, max=d'arco)

# Inner loop voltage -
kpv_range = collect(0.5:0.03:0.6) # D'Arco range 
kiv_range = collect(400:100:800) # D'Arco range 

# Inner loop current gains 
kpc_range = collect(0.74:0.15:1.27) #? min=markovic max=d'arco
kic_range = collect(1.19:4.0:14.3) #? min=markovic, max=d'arco

load_scale_range = collect(0.5:0.3:1.5)
inverter_share_range = collect(0.0:0.2:1.0)
gfm_share_range = collect(0.0:0.2:1.0)

col_names = ["load_scale", "kq", "kpv", "kiv", "kpc", "kic", "inv_share", "gfm_share"]
df = DataFrame([name => [] for name in col_names])

for kq = gfm_kq_range;
    for kpv = kpv_range;
        for kiv = kiv_range;
            for kpc = kpc_range;
                for kic = kic_range;
                    for load_scale = load_scale_range;
                        for inv_share = inverter_share_range;
                            for gfm_share = gfm_share_range;
                                push!(df, [load_scale, kq, kpv, kiv, kpc, kic, inv_share, gfm_share])
                            end
                        end
                    end
                end
            end
        end
    end
end

CSV.write("nrel_parameters_kpv_kiv_adjusted.csv", df)