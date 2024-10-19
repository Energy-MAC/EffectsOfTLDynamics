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

gss = load_serde_data("data/gab_tests")
df = gss.df

df.Case = df."injector at {Bus1}".*" ".*df."injector at {Bus 2}".*" ".*df."injector at {Bus 3}"

filter!(row -> row.sim.results != nothing, df)

# filter(row -> row."Line Model" == "statpi", df)

# filter(row -> row."Line scale" >= 2.0 , df)


add_result!(gss,
    ["Bus 1 Injector Current", "Bus 2 Injector Current", "Bus 3 Injector Current"],
    PSE.get_injector_currents,
)
# tell it to record the timestamps
add_result!(gss, "Time", PSE.get_time)
# add_result!(gss, "initial_sm", PSE.get_sm)
# add_result!(gss, "final_sm", PSE.small_signal_tripped)
# add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)

add_result!(gss, "Eigs", PSE.get_eigenvalues)
sort!(df, [:Case, :"Line scale"])

p = makeplots(
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
    xaxis_home_range = (min=0.48, max=0.55)  
)

savehtmlplot(p, "first.html")

p = makeplots(df, 
    supertitle = "Transient results",
    x = "Time",
    y = "Bus 2 Injector Current",
    x_title = "Time [s]",
    y_title = "Current at Bus 2",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    # legendgroup = "Line Model",
    color = "Line Model",
    scattermode = "lines",
    yaxis_home_range = (min=0.0, max=2.0),
    xaxis_home_range = (min=0.48, max=0.55)  
)

savehtmlplot(p, "second.html")

p = makeplots(df, 
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
    xaxis_home_range = (min=0.48, max=0.55)  
)

savehtmlplot(p, "third.html")