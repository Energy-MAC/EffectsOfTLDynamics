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

pathname = "../../../../data/gabrielecr/"*line_params.alg_line_name*" line trip/"*current_time_string
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
sort!(df, [:"Line scale", :"Load scale", :"Line Model"])

# selected_df = df[:, [:"Line Model", :"Line scale", :"Load scale", :"Real", :"Imag", :"Max real part"]]
# using CSV
# CSV.write("output.csv", selected_df)

# filter!(row -> row.Case == "SM SM GFL", df)

# df."real" = map(x->x[1], df."Real")
# df."imag" = map(x->x[1], df."Imag")
# save_serde_data(gss, "data/gab_tests")
# sort!(df, [:"Line scale", :"Load scale", :"Line Model"])

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
    xaxis_home_range = (min=0.48, max=2.0)  
)

pathname2 = pathname*"/plots"
savehtmlplot(p, pathname2)

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
    xaxis_home_range = (min=0.48, max=2.0)  
)

savehtmlplot(p, "i_bus2_2.html")

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
    xaxis_home_range = (min=0.48, max=2.0)  
)

savehtmlplot(p, "i_bus3_2.html")

p = makeplots(df, 
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

savehtmlplot(p, "eigs_2.html")

p = makeplots(
    df, 
    supertitle = "Transient results",
    x = "Time",
    y = "Voltage at Bus 5",
    x_title = "Time [s]",
    y_title = "Voltage at Bus 4",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    legendgroup = "Line Model",
    color = "Line Model",
    scattermode = "lines",
    yaxis_home_range = (min=0.0, max=2.0),
    xaxis_home_range = (min=0.48, max=2.0)  
)

savehtmlplot(p, "v_bus5.html")