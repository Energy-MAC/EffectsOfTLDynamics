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


line_scales = collect(1.0:4.0:5.0)
load_scales = collect(0.5:0.5:1.0)

vars_to_measure = ["inverter_currents.csv"]

plots = []
plt = []

line_scales = [1.0]
load_scales = [0.5]

for var_to_measure in vars_to_measure
    for line_scale in line_scales
        for load_scale in load_scales
            
            partial_path = main_path*"/$(line_scale)_$(load_scale)"
            folder_path = partial_path*"/statpi"
            df_alg = CSV.read(partial_path*"/statpi/"*var_to_measure, DataFrame);
            df_dyn = CSV.read(partial_path*"/dynpi/"*var_to_measure, DataFrame); 
            df_ms = CSV.read(partial_path*"/MSSB/"*var_to_measure, DataFrame);
            df_msmb = CSV.read(partial_path*"/MSMB/"*var_to_measure, DataFrame);

            alg_idx = findall(df_alg."Time" .< 0.275)[end]
            dyn_idx = findall(df_dyn."Time" .< 0.275)[end]
            ms_idx = findall(df_ms."Time" .< 0.275)[end]
            msmb_idx = findall(df_msmb."Time" .< 0.275)[end]

            sol_alg = (df_alg."Time"[1:alg_idx], df_alg."generator-102-1_i_mag"[1:alg_idx])
            sol_dyn = (df_dyn."Time"[1:dyn_idx], df_dyn."generator-102-1_i_mag"[1:dyn_idx])
            sol_ms = (df_ms."Time"[1:ms_idx], df_ms."generator-102-1_i_mag"[1:ms_idx])
            sol_msmb = (df_msmb."Time"[1:msmb_idx], df_msmb."generator-102-1_i_mag"[1:msmb_idx])
        
            plt = plot(sol_alg, label = L"\mathrm{statpi}", legend = :bottomright)
            plot!(plt, sol_dyn, label = L"\mathrm{dynpi}")
            plot!(plt, sol_ms, label = L"\mathrm{MSSB}")
            plot!(plt, sol_msmb, label = L"\mathrm{MSMB}")
            plot!(plt, legend = true)    
            plot!(xlabel = L"$ \mathrm{Time} \quad [s]$", title = "")  
            plot!(ylabel = L"$||I_f|| \quad \mathrm{[\ p.u.]}$")
            plot!(framestyle=:box)
            plot!(size=(1000,800))
            plot!(xguidefontsize=25, yguidefontsize=25, tickfontsize=16,legendfontsize=16)
            
            savefig(partial_path*"/"*var_to_measure[1:end-4]*"_mags.svg")
            plot!(xlims=(0.249, 0.275))
            savefig(partial_path*"/"*var_to_measure[1:end-4]*"_mags_100ms_zoom.svg")
            # push!(plots, plt)
        end 
    end
end