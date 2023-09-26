cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics

using CSV
using DataFrames
using Dates
using LaTeXStrings
using Plots

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

rn = "2023-09-24T23:14:16.337" #string date
model = "9bus_slackless"
main_path = "../results/"*model*"/"*rn

df_alg = nothing
df_dyn = nothing
df_ms = nothing
df_msmb = nothing

sol_alg = nothing
sol_dyn = nothing
sol_ms = nothing
sol_msmb = nothing


line_scales = collect(1.0:1:3.0)
load_scales = collect(0.5:0.5:2.0)

line_scales = collect(3.0)
load_scales = collect(0.5:0.5:1)

vars_to_measure = ["inverter_currents.csv"]

plots = []
plt = []

for var_to_measure in vars_to_measure
    for line_scale in line_scales
        for load_scale in load_scales
            
            partial_path = main_path*"/$(line_scale)_$(load_scale)"
            folder_path = partial_path*"/statpi"
            
            df_alg = CSV.read(partial_path*"/statpi/"*var_to_measure, DataFrame);
            ir_alg = df_alg."generator-1-1_ir_cnv"
            ii_alg = df_alg."generator-1-1_ii_cnv"
            i_alg_mag = sqrt.(ir_alg.^2 + ii_alg.^2)
            df_alg[!,"generator-1-1_i_mag"] = i_alg_mag
            CSV.write(partial_path*"/statpi/"*var_to_measure, df_alg)
            
            df_dyn = CSV.read(partial_path*"/dynpi/"*var_to_measure, DataFrame); 
            ir_dyn = df_dyn."generator-1-1_ir_cnv"
            ii_dyn = df_dyn."generator-1-1_ii_cnv"
            i_dyn_mag = sqrt.(ir_dyn.^2 + ii_dyn.^2)
            df_dyn[!,"generator-1-1_i_mag"] = i_dyn_mag
            CSV.write(partial_path*"/dynpi/"*var_to_measure, df_dyn)

            df_ms = CSV.read(partial_path*"/MSSB/"*var_to_measure, DataFrame);
            ir_ms = df_ms."generator-1-1_ir_cnv"
            ii_ms = df_ms."generator-1-1_ii_cnv"
            i_ms_mag = sqrt.(ir_ms.^2 + ii_ms.^2)
            df_ms[!,"generator-1-1_i_mag"] = i_ms_mag
            CSV.write(partial_path*"/MSSB/"*var_to_measure, df_ms)

            df_msmb = CSV.read(partial_path*"/MSMB/"*var_to_measure, DataFrame);
            ir_msmb = df_msmb."generator-1-1_ir_cnv"
            ii_msmb = df_msmb."generator-1-1_ii_cnv"
            i_msmb_mag = sqrt.(ir_msmb.^2 + ii_msmb.^2)
            df_msmb[!,"generator-1-1_i_mag"] = i_msmb_mag
            CSV.write(partial_path*"/MSMB/"*var_to_measure, df_msmb)

        end 
    end
end
