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
function add_results_gss!()
    add_result!(gss,
        ["Bus 1 Injector Current", "Bus 2 Injector Current", "Bus 3 Injector Current"],
        PSE.get_injector_currents,
    )
    add_result!(gss, "Time", PSE.get_time)
    add_result!(gss, "Eigs", PSE.get_eigenvalues)
end
add_results_gss!()

df = gss.df
df.Case = df."injector at {Bus1}".*" ".*df."injector at {Bus 2}".*" ".*df."injector at {Bus 3}"
df."Line scale" = (x->x.line_scale).(df."Line Params")
df."Eigs"
df."Real" = (x->real(x)).(df."Eigs")
df."real" = map(x->x[1], df."Real")
df."Imag" = (x->imag(x)).(df."Eigs")
df."imag" = map(x->x[1], df."Imag")
# save_serde_data(gss, "data/gab_tests")

# filter!(row -> row.sim.results != no(hing, df)
# filter(row -> row."Line Model" == "statpi", df)
# filter(row -> row."Line scale" >= 2.0 , df)
# add_result!(gss,
#     ["Bus 1 Injector Current", "Bus 2 Injector Current", "Bus 3 Injector Current"],
#     PSE.get_injector_currents,
# )
# tell it to record the timestamps
# add_result!(gss, "Time", PSE.get_time)
# add_result!(gss, "initial_sm", PSE.get_sm)
# add_result!(gss, "final_sm", PSE.small_signal_tripped)
# add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)

sort!(df, [:Case, :"Line scale", :"Line Model"])

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

savehtmlplot(p, "bus1.html")

p = makeplots(df, 
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
    xaxis_home_range = (min=0.48, max=0.55)  
)

savehtmlplot(p, "bus2.html")

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

savehtmlplot(p, "bus3.html")

p = makeplots(df, 
    supertitle = "Small signal results",
    x = "real",
    y = "imag",
    x_title = "Real",
    y_title = "Imag",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    legendgroup = "Line Model",
    color = "Line Model",
    scattermode = "markers",
    # yaxis_home_range = (min=0.0, max=2.0),
    # xaxis_home_range = (min=0.48, max=0.55)  
)

savehtmlplot(p, "eigs.html")