cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "../raw_data/WSCC_9bus.raw"))

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

#Define machine
# Create the machine
machine_oneDoneQ() = OneDOneQMachine(
    0.0, #R
    1.3125, #Xd
    1.2578, #Xq
    0.1813, #Xd_p
    0.25, #Xq_p
    5.89, #Td0_p
    0.6, #Tq0_p
)

# Shaft
shaft_no_damping() = SingleMass(
    3.01, #H (M = 6.02 -> H = M/2)
    0.0, #D
)

# AVR: Type I: Resembles a DC1 AVR
avr_type1() = AVRTypeI(
    20.0, #Ka - Gain
    0.01, #Ke
    0.063, #Kf
    0.2, #Ta
    0.314, #Te
    0.35, #Tf
    0.001, #Tr
    (min = -5.0, max = 5.0),
    0.0039, #Ae - 1st ceiling coefficient
    1.555, #Be - 2nd ceiling coefficient
)

#No TG
tg_none() = TGFixed(1.0) #efficiency

#No PSS
pss_none() = PSSFixed(0.0) #Vs

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
    #Find the generator at bus 102
    if get_number(get_bus(g)) == 2
        #Create the dynamic generator
        case_gen = DynamicGenerator(
            get_name(g),
            1.0, # ω_ref,
            machine_oneDoneQ(), #machine
            shaft_no_damping(), #shaft
            avr_type1(), #avr
            tg_none(), #tg
            pss_none(), #pss
        )
        #Attach the dynamic generator to the system by
        # specifying the dynamic and static components
        add_component!(sys, case_gen, g)
        #Find the generator at bus 103
    elseif get_number(get_bus(g)) == 3
        #Create the dynamic inverter
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
        #Attach the dynamic inverter to the system
        add_component!(sys, case_inv, g)
    end
end

to_json(sys, joinpath(pwd(), "../json_data/9bus.json"), force = true)