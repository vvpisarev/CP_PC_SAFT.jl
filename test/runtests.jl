using Test
using CubicEoS
using CP_PC_SAFT

mixture = CubicEoS.load(CPPCSAFTMixture; names=("methane", "n-pentane"))
volume = 1e-6
nmol = [0.6, 0.4] * 7000 * volume
RT = CubicEoS.GAS_CONSTANT_SI * 500

@test pressure(mixture, nmol, volume, RT) ≈ 2.8631703619799692e7
