cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
using EffectsOfTLDynamics
using ControlSystems

using CSV
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

### Extract line data from files for M=1 and M=3
M = 1
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

M = 3
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

mssb_r = real(z_km_1)[1]
mssb_l = imag(z_km_1[1])/(2*pi*60)
msmb_rs = real(z_km_3)
msmb_ls = imag(z_km_3)/(2*pi*60)

# Make Z(s) plots
mssb_tf = tf([1], [mssb_l, mssb_r])

msmb_tf = 0;
for i = [1,2,3];
    msmb_tf = msmb_tf + tf(1, [msmb_ls[i], msmb_rs[i]]);
end

bodeplot(msmb_tf, label="MSMB")
bodeplot!(mssb_tf, label="MSSB")
xlims!(10^0, 10^3)

bodeplot(mssb_tf-msmb_tf, label="SB minus MB")
