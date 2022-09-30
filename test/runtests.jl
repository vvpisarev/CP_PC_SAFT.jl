using Test
using CubicEoS
using CP_PC_SAFT

@testset "CP_PC_SAFT.jl" begin

c1c5 = CubicEoS.load(CPPCSAFTMixture; names=("methane", "n-pentane"))
volume = 1e-6
nmol = [0.6, 0.4] * 7000 * volume
RT500 = CubicEoS.GAS_CONSTANT_SI * 500

@testset "basic_thermo" begin
    @test pressure(c1c5, nmol, volume, RT500) ≈ 2.8631703619799692e7
end

end # @testset
