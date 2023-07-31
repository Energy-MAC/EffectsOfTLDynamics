using PowerSystems
using PowerSimulationsDynamics

include("experiment_structs.jl")

default_bic_params = BICParam("line-102-103", 1.0)
default_gen_trip_params = GenTripParam(DynamicGenerator, "generator-102-1")
default_crc_params = CRCParam(DynamicGenerator, "generator-102-1", :V_ref, 0.95)
default_lc_params = LCParam(ElectricLoad, "load-103-1", :P_ref, 1.08)
default_lt_params = LTParam(ElectricLoad, "load-103-1")
default_sbvc_params = SBVCParam(:V_ref, 1.048)

export default_bic_params
export default_gen_trip_params
export default_crc_params
export default_lc_params
export default_lt_params
export default_sbvc_params
