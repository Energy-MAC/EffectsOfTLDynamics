cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using ZIPE_loads
using InfrastructureSystems
using PlotlyJS
using DataFrames
using TLmodels
using DataFramesMeta
using LaTeXStrings
using Colors

using PowerSystemsExperiments
const PSE = PowerSystemsExperiments

# pathname = "../../../../data/gabrielecr/Bus 5-Bus 4-i_1/2024-11-27T18:32:38.427"
pathname = "../../../../data/gabrielecr/Bus 7-Bus 5-i_1/2024-11-27T17:20:32.588"
pathname = "../../../../data/gabrielecr/Bus 7-Bus 5-i_1/2024-12-06T01:04:59.235"

pathname = "../../../../data/gabrielecr/BUS 1-BUS 2-i_1_static/2024-12-06T04:01:45.783"

pathname2 = pathname*"/plots"
mkdir(pathname2)

gss = load_serde_data(pathname)

function add_results_gss!()
    add_result!(gss,
        ["Bus 1 Injector Current", "Bus 2 Injector Current"],
        PSE.get_injector_currents,
    )
    add_result!(gss, "Time", PSE.get_time)
    add_result!(gss, ["Voltage at $busname" for busname in get_name.(get_components(ACBus, gss.base))], PSE.get_bus_voltages)
    add_result!(gss, "Eigs", PSE.get_eigenvalues)
end
add_results_gss!()

df = gss.df
df.Case = df."injector at {BUS 1}".*" ".*df."injector at {BUS 2}"
df."Line scale" = (x->x.line_scale).(df."Line Params")
df."Eigs"
df."Real" = (x->real(x)[1]).(df."Eigs")
df."Imag" = (x->imag(x)[1]).(df."Eigs")
df."Max real part" = maximum.(df."Real")
sort!(df, [:"Line scale", :"Load scale", :"Line Model", :"Case"])