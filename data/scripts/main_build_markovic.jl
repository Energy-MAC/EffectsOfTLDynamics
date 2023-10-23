cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using InfrastructureSystems

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "../raw_data/OMIB.raw"))

# Define machine
# Create the machine
AF_machine() = AndersonFouadMachine(
    0.0041, # R::Float64
    1.25, # Xd::Float64
    1.22, # Xq::Float64
    0.232, # Xd_p::Float64
    0.715, # Xq_p::Float64
    0.12, # Xd_pp::Float64
    0.12, # Xq_pp::Float64
    4.75, # Td0_p::Float64
    1.5, # Tq0_p::Float64
    0.06, # Td0_pp::Float64
    0.21# Tq0_pp::Float64
)

# Shaft
shaft_no_damping() = SingleMass(
    6.5, # H (M = 6.02 -> H = M/2)
    1.0, #0  #D
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
tg_none() = TGFixed(1.0)

# TGTypeI(
#     0.02,
#     0.1,
#     0.45,
#     0.0,
#     12.0,
#     50.0,
#     (min = 0.3, max = 1.2)
# ) #efficiency

#No PSS
pss_none() = STAB1(
    5, # KT
    2, #T
    8.3, #T1T3
    0.03, #T3
    10, #T2T4
    0.015, #T4
    0.2
)

#Define converter as an AverageConverter
converter_high_power() = AverageConverter(
    rated_voltage = 138.0, 
    rated_current = 100.0
    )

function droop_outer_control()
    function active_droop()
        return PSY.ActivePowerDroop(; Rp = 0.02, ωz = 2 * pi * 5)
    end
    function reactive_droop()
        return ReactivePowerDroop(; kq = 0.001, ωf = 2 * pi * 5)
    end
    return OuterControl(active_droop(), reactive_droop())
end


#Define an Inner Control as a Voltage+Current Controler with Virtual Impedance:
GFM_inner_control() = VoltageModeControl(
    kpv = 0.52,     #Voltage controller proportional gain
    kiv = 1.16,    #Voltage controller integral gain
    kffv = 1.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
    rv = 0.001,       #Virtual resistance in pu
    lv = 0.2,       #Virtual inductance in pu
    kpc = 0.74,     #Current controller proportional gain
    kic = 1.19,     #Current controller integral gain
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
filt() = LCLFilter(lf = 0.08, rf = 0.03, cf = 0.074, lg = 0.2, rg = 0.01)

# build static component for gen1
gen2 = get_component(Generator, sys, "generator-102-1")
gen1 = deepcopy(gen2);
gen1.name = "generator-101-1";
bus1 = get_component(Bus, sys, "BUS 1")
gen1.bus = bus1
gen1.time_series_container = InfrastructureSystems.TimeSeriesContainer(); 
# Add gen1 to sys
add_component!(sys, gen1)

# Droop GFM v Machine
for g in get_components(Generator, sys)
    if get_number(get_bus(g)) == 101;
        case_gen = DynamicGenerator(
            get_name(g),
            1.0, # ω_ref,
            AF_machine(), #machine
            shaft_no_damping(), #shaft
            avr_type1(), #avr
            tg_none(), #tg
            pss_none(), #pss
        )
        #case_gen.bus.BusTypes = BusTypes.REF
        #Attach the dynamic generator to the system by
        # specifying the dynamic and static components
        add_component!(sys, case_gen, g)
    elseif get_number(get_bus(g)) == 102;
        case_inv = DynamicInverter(
                    get_name(g),
                    1.0, # ω_ref,
                    converter_high_power(), #converter
                    droop_outer_control(), #outer control
                    GFM_inner_control(), #inner control voltage source
                    dc_source_lv(), #dc source
                    pll(), #pll
                    filt(), #filter
                )
        add_component!(sys, case_inv, g)
    end
end

# Save as json 
to_json(sys, joinpath(pwd(), "../json_data/gorno.json"), force = true)
