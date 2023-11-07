using Test
using CubicEoS
using CP_PC_SAFT

@testset "CP_PC_SAFT.jl" begin

mixture = CubicEoS.load(CPPCSAFTMixture; names=("methane", "n-pentane"))
volume = 1e-6
nmol = [0.6, 0.4] * 7000 * volume
RT = CubicEoS.GAS_CONSTANT_SI * 500

@testset "basic_thermo" begin
    @test pressure(mixture, nmol, volume, RT) ≈ 2.8631703619799692e7
end

@testset "vt_stability" begin
    @testset "Stable one phase state" begin
        issuccess, isstable, trials = vt_stability(mixture, nmol, volume, RT;
            tol=1e-6,
            maxiter=1000,
        )
        @test issuccess && isstable
    end

    @testset "Unstable one phase state" begin
        RTUnstable = CubicEoS.GAS_CONSTANT_SI * 300
        issuccess, isstable, trials = vt_stability(mixture, nmol, volume, RTUnstable;
            tol=1e-6,
            maxiter=1000,
        )
        @test issuccess && !isstable
    end
end

@testset "vt_split" begin
    RTUnstable = CubicEoS.GAS_CONSTANT_SI * 300
    _, _, trials = vt_stability(mixture, nmol, volume, RTUnstable;
        tol=1e-6,
        maxiter=1000,
    )
    conc = CubicEoS.concentrationwithlowesttpd(trials)
    result = vt_split(mixture, nmol, volume, RTUnstable, conc, CubicEoS.VTSplitIdealIdentityState)
    @test converged(result)
end

end # @testset
