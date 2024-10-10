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

get_inv_d_fluxlink = load_serde_data("data/gab_tests/results0.jls")

df."Case" = [df."injector at {Bus1}".*" ".*df."injector at {Bus 2}".*" ".*df."injector at {Bus 3}"]

p = makeplots(df, 
    x = "Time",
    y = "Bus 1 Injector Current",
    rows = "Line scale",
    cols = "Case",
    slider = "Load scale",
    legendgroup = "Line Model",
    color = "Line Model",
    scattermode = "lines",
    x_title = "Time [s]",
    y_title = "Current at Bus 1",
    supertitle = "Transient results"
)
