using PowerSystems
using PowerSimulationsDynamics
using Sundials
using Plots
using PowerNetworkMatrices
using SparseArrays
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "test_sys.json"));
sys_iter = deepcopy(sys)

tspan = (0.0, 30.0)

N = 20
for l in get_components(Line, sys_iter)
    bus_from = l.arc.from
    bus_to = l.arc.to
    # Create a bunch of Bus
    start_bus = bus_from
    for b_ix in 1:N - 1
        println(b_ix)
        bus_to_create = Bus(
            number = 1000000000 + 100000*bus_from.number + 100*bus_to.number + b_ix,
            name = bus_from.name * "-" * bus_to.name * "-internal-bus_" * string(b_ix),
            bustype = BusTypes.PQ,
            angle = 0.0,
            magnitude = 1.0,
            voltage_limits = (min = 0.95, max = 1.05),
            base_voltage = bus_from.base_voltage,
        )
        add_component!(sys, bus_to_create)
        end_bus = bus_to_create
        line_to_create = Line(
            name = l.name * "_segment_" * string(b_ix),
            available = true,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = Arc(from = start_bus, to = end_bus),
            r = l.r/N,
            x = l.x/N,9
            b = (from = l.b.from * N, to = l.b.to * N),
            rate = l.rate,
            angle_limits = (min = -pi/2, max = pi/2),
        )
        add_component!(sys, line_to_create)
        println(get_name(line_to_create))
        start_bus = end_bus
    end
    line_to_create = Line(
            name = l.name * "_segment_" * string(N),
            available = true,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = Arc(from = start_bus, to = bus_to),
            r = l.r/N,
            x = l.x/N,
            b = (from = l.b.from * N, to = l.b.to * N),
            rate = l.rate,
            angle_limits = (min = -pi/2, max = pi/2),
    )
    add_component!(sys, line_to_create)
    println(get_name(line_to_create))
    remove_component!(sys, l)
    remove_component!(sys, l.arc)
end
show_components(sys,Line)
show_components(sys, Bus)

dyn_sys = deepcopy(sys)

using NLsolve

for l in get_components(Line, sys)
    print(l.name)
end

function f(x,p)
    γ = sqrt(x[1]*x[2])
    return [z - x[1]*l*(sinh(γ*l)/(γ*l)), y - x[2]*l/2*(tanh(γ*l)/(γ*l))]
end

line = get_component(Line, sys, "BUS 1-BUS 2-i_1")
z = line.r + im*line.x
y = im*(line.b.from + line.b.to)
l = 100 #km, defined outside
sol = nlsolve(f, [0.01 + im*0.12, 0 + im*0.06])
if sol.x_converged == false
    error("Impedance calc error")
end
z_km, y_km = sol.zero

using NonlinearSolve, StaticArrays
u0 = @SVector[0.01 + im*0.5, 0 + im*3]
p = 0
probN = NonlinearProblem(f, u0, p)
solver = solve(probN, NewtonRaphson(), reltol = 1e-9)

z_line = z_km * l
y_line = y_km * l

line_to_create = Line(
    name = l.name * "_segment_" * string(b_ix),
    available = true,
    active_power_flow = 0.0,
    reactive_power_flow = 0.0,
    arc = Arc(from = start_bus, to = end_bus),
    r = real(z_line)/N,
    x = imag(z_line)/N,
    b = (from = imag(y_line) * N, to = imag(y_line) * N),
    rate = l.rate,
    angle_limits = l.angle_limits,
)
#230kV line @ 60 Hz
Z_c = 380 
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
z_km = r_km + im*x_km # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
y_km = g_km + im*b_km
z_km_pu = z_km/Z_c
y_km_pu = y_km*Z_c
l = 100 #km
γ = sqrt(z_km*y_km)
z_ll = z_km_pu*l*(sinh(γ*l)/(γ*l))
y_ll = y_km_pu/2*l*(tanh(γ*l)/(γ*l))