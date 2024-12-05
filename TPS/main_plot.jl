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

pathname2 = pathname*"/plots"
mkdir(pathname2)

gss = load_serde_data(pathname)

function add_results_gss!()
    add_result!(gss,
        ["Bus 1 Injector Current", "Bus 2 Injector Current", "Bus 3 Injector Current"],
        PSE.get_injector_currents,
    )
    add_result!(gss, "Time", PSE.get_time)
    add_result!(gss, ["Voltage at $busname" for busname in get_name.(get_components(ACBus, gss.base))], PSE.get_bus_voltages)
    add_result!(gss, "Eigs", PSE.get_eigenvalues)
end
add_results_gss!()

df = gss.df
df.Case = df."injector at {Bus1}".*" ".*df."injector at {Bus 2}".*" ".*df."injector at {Bus 3}"
df."Line scale" = (x->x.line_scale).(df."Line Params")
df."Eigs"
df."Real" = (x->real(x)[1]).(df."Eigs")
df."Imag" = (x->imag(x)[1]).(df."Eigs")
df."Max real part" = maximum.(df."Real")
sort!(df, [:"Line scale", :"Load scale", :"Line Model", :"Case"])

# selected_df = df[:, [:"Line Model", :"Line scale", :"Load scale", :"Real", :"Imag", :"Max real part"]]
# using CSV
# CSV.write("output.csv", selected_df)

# filter!(row -> row.Case == "SM SM GFL", df)

# df."real" = map(x->x[1], df."Real")
# df."imag" = map(x->x[1], df."Imag")
# save_serde_data(gss, "data/gab_tests")
# sort!(df, [:"Line scale", :"Load scale", :"Line Model"])

p1 = makeplots(
    df, 
    supertitle = "Transient results",
    x = "Time",
    y = "Bus 1 Injector Current",
    x_title = "Time [s]",
    y_title = "Current at Bus 1",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    legendgroup = "Line Model",
    color = "Line Model",
    scattermode = "lines",
    yaxis_home_range = (min=0.0, max=2.0),
    xaxis_home_range = (min=0.48, max=2.0)  
)

savehtmlplot(p1, pathname2*"/i_bus1.html")

p2 = makeplots(df, 
    supertitle = "Transient results",
    x = "Time",
    y = "Bus 2 Injector Current",
    x_title = "Time [s]",
    y_title = "Current at Bus 2",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    legendgroup = "Line Model",
    color = "Line Model",
    scattermode = "lines",
    yaxis_home_range = (min=0.0, max=2.0),
    xaxis_home_range = (min=0.48, max=2.0)  
)

savehtmlplot(p2, pathname2*"/i_bus2.html")

p3 = makeplots(df, 
    supertitle = "Transient results",
    x = "Time",
    y = "Bus 3 Injector Current",
    x_title = "Time [s]",
    y_title = "Current at Bus 3",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    legendgroup = "Line Model",
    color = "Line Model",
    scattermode = "lines",
    yaxis_home_range = (min=0.0, max=2.0),
    xaxis_home_range = (min=0.48, max=2.0)  
)

savehtmlplot(p3, pathname2*"/i_bus3.html")

pe = makeplots(df, 
    supertitle = "Small signal results",
    x = "Real",
    y = "Imag",
    x_title = "Real",
    y_title = "Imag",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    legendgroup = "Line Model",
    color = "Line Model",
    #scattermode = "markers",
)

savehtmlplot(pe, pathname2*"/eigs.html")