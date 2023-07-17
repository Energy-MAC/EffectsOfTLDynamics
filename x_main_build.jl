using PowerSystems
using PowerSimulationsDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "OMIB_new.raw"))

slack_bus = [b for b in get_components(Bus, sys) if get_bustype(b) == BusTypes.REF][1]

inf_source = Source(
           name = "InfBus", #name
           available = true, #availability
           active_power = 0.0,
           reactive_power = 0.0,
           bus = slack_bus, #bus
           R_th = 0.0, #Rth
           X_th = 5e-6, #Xth
       )

add_component!(sys, inf_source)

#Define converter as an AverageConverter
converter_high_power() = AverageConverter(rated_voltage = 138.0, rated_current = 100.0)

#Define Outer Control as a composition of Virtual Inertia + Reactive Power Droop
outer_control() = OuterControl(
    VirtualInertia(Ta = 2.0, kd = 400.0, kω = 20.0),
    ReactivePowerDroop(kq = 0.2, ωf = 1000.0),
)

#Define an Inner Control as a Voltage+Current Controler with Virtual Impedance:
inner_control() = VoltageModeControl(
    kpv = 0.59,     #Voltage controller proportional gain
    kiv = 736.0,    #Voltage controller integral gain
    kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
    rv = 0.0,       #Virtual resistance in pu
    lv = 0.2,       #Virtual inductance in pu
    kpc = 1.27,     #Current controller proportional gain
    kic = 14.3,     #Current controller integral gain
    kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
    ωad = 50.0,     #Active damping low pass filter cut-off frequency
    kad = 0.2,      #Active damping gain
)

#Define DC Source as a FixedSource:
dc_source_lv() = FixedDCSource(voltage = 600.0)

#Define a Frequency Estimator as a PLL based on Vikram Kaura and Vladimir Blaskoc 1997 paper:
pll() = KauraPLL(
    ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
    kp_pll = 0.084,  #PLL proportional gain
    ki_pll = 4.69,   #PLL integral gain
)

#Define an LCL filter:
filt() = LCLFilter(lf = 0.08, rf = 0.003, cf = 0.074, lg = 0.2, rg = 0.01)
#filt() = LCFilter(lf = 0.08, rf = 0.003, cf = 0.074)

for g in get_components(Generator, sys)
    case_inv = DynamicInverter(
        "generator-102-1",
        1.0, # ω_ref,
        converter_high_power(), #converter
        outer_control(), #outer control
        inner_control(), #inner control voltage source
        dc_source_lv(), #dc source
        pll(), #pll
        filt(), #filter
    )

#Attach the dynamic inverter to the system
    add_component!(sys, case_inv, g)
end

to_json(sys, joinpath(pwd(), "test_sys.json"), force = true)