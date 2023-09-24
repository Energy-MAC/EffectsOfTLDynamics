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

rn = "2023-09-24T20:12:13.348" #string date
model = "inv_v_machine"
main_path = "../results/"*model*"/"*rn

df_alg = nothing
df_dyn = nothing
df_ms = nothing
df_msmb = nothing

sol_alg = nothing
sol_dyn = nothing
sol_ms = nothing
sol_msmb = nothing


line_scales = collect(1.0:2.0:7.0)
load_scales = collect(0.5:0.5:3.0)

line_scales = collect(7.0)
load_scales = collect(0.5:0.5:1.5)

vars_to_measure = ["inverter_currents.csv"]

plots = []
plt = []

for var_to_measure in vars_to_measure
    for line_scale in line_scales
        for load_scale in load_scales
            
            partial_path = main_path*"/$(line_scale)_$(load_scale)"
            folder_path = partial_path*"/statpi"
            df_alg = CSV.read(partial_path*"/statpi/"*var_to_measure, DataFrame);
            df_dyn = CSV.read(partial_path*"/dynpi/"*var_to_measure, DataFrame); 
            df_ms = CSV.read(partial_path*"/MSSB/"*var_to_measure, DataFrame);
            df_msmb = CSV.read(partial_path*"/MSMB/"*var_to_measure, DataFrame);
            sol_alg = (df_alg."Time", df_alg."generator-102-1_ir_cnv")
            sol_dyn = (df_dyn."Time", df_dyn."generator-102-1_ir_cnv")
            sol_ms = (df_ms."Time", df_ms."generator-102-1_ir_cnv")
            sol_msmb = (df_msmb."Time", df_msmb."generator-102-1_ir_cnv")
        
            plt = plot(sol_alg, label = L"\mathrm{statpi}")
            plot!(plt, sol_dyn, label = L"\mathrm{dynpi}")
            plot!(plt, sol_ms, label = L"\mathrm{MSSB}")
            plot!(plt, sol_msmb, label = L"\mathrm{MSMB}")
            plot!(plt, legend = true)    
            plot!(xlabel = L"$ \mathrm{Time} \quad [s]$", title = "")  
            plot!(ylabel = L"$||V_2|| \quad \mathrm{[\ p.u.]}$")
            plot!(legend_title = L"$%$load_scale \ %$line_scale$")
            
            savefig(partial_path*"/"*var_to_measure[1:end-4]*".svg")

            plot!(xlims=(0.249, 0.26))
            savefig(partial_path*"/"*var_to_measure[1:end-4]*"_100ms_zoom.svg")
            # push!(plots, plt)
        end 
    end
end
