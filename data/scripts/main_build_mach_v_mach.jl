cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using InfrastructureSystems

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "../raw_data/OMIB.raw"))

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

# build static component for gen1
gen2 = get_component(Generator, sys, "generator-102-1")
gen1 = deepcopy(gen2);
gen1.name = "generator-101-1";
bus1 = get_component(Bus, sys, "BUS 1")
gen1.bus = bus1
gen1.time_series_container = InfrastructureSystems.TimeSeriesContainer(); 
# Add gen1 to sys
add_component!(sys, gen1)

for g in get_components(Generator, sys)
    if get_number(get_bus(g)) == 101;
        case_gen = DynamicGenerator(
            get_name(g),
            1.0, # ω_ref,
            machine_oneDoneQ(), #machine
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
        case_gen = DynamicGenerator(
            get_name(g),
            1.0, # ω_ref,
            machine_oneDoneQ(), #machine
            shaft_no_damping(), #shaft
            avr_type1(), #avr
            tg_none(), #tg
            pss_none(), #pss
        )
        add_component!(sys, case_gen, g)
    end

end

# Save as json 
to_json(sys, joinpath(pwd(), "../json_data/mach_v_mach.json"), force = true)
