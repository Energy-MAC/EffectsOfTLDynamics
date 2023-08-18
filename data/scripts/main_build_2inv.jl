cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "../raw_data/OMIB.raw"))

# Copy the bus 2 generator
gen2 = get_component(Generator, sys, "generator-102-1")
gen1 = deepcopy(gen2);
gen1.name = "generator-101-1";
bus1 = get_component(Bus, sys, "BUS 1")
gen1.bus = bus1
gen1.time_series_container = InfrastructureSystems.TimeSeriesContainer(); 
# Add gen1 to sys
add_component!(sys, gen1)


#Define converter as an AverageConverter
converter_high_power() = AverageConverter(
    rated_voltage = 138.0, 
    rated_current = 100.0
    )

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
dc_source_lv() = FixedDCSource(
    voltage = 600.0
    )

#Define a Frequency Estimator as a PLL based on Vikram Kaura and Vladimir Blaskoc 1997 paper:
pll() = KauraPLL(
    ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
    kp_pll = 0.084,  #PLL proportional gain
    ki_pll = 4.69,   #PLL integral gain
)

#Define an LCL filter:
filt() = LCLFilter(lf = 0.08, rf = 0.003, cf = 0.074, lg = 0.2, rg = 0.01)


for g in get_components(Generator, sys)
    if get_number(get_bus(g)) == 101;
        case_inv = DynamicInverter(
                    get_name(g),
                    1.0, # ω_ref,
                    converter_high_power(), #converter
                    outer_control(), #outer control
                    inner_control(), #inner control voltage source
                    dc_source_lv(), #dc source
                    pll(), #pll
                    filt(), #filter
                )
        add_component!(sys, case_inv, g)
    elseif get_number(get_bus(g)) == 102;
        case_inv = DynamicInverter(
                    get_name(g),
                    1.0, # ω_ref,
                    converter_high_power(), #converter
                    outer_control(), #outer control
                    inner_control(), #inner control voltage source
                    dc_source_lv(), #dc source
                    pll(), #pll
                    filt(), #filter
                )
        add_component!(sys, case_inv, g)
    end

end

# Save as json 
to_json(sys, joinpath(pwd(), "../json_data/twobus_2inv.json"), force = true)
