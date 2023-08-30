using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test


lode = lodeproblem()
pref = exact_solution(podeproblem())


@testset "$(rpad("Hamilton-Pontryagin integrators",80))" begin

    ϕ(v, q̄, q, a, Δt) = v .= (q .- q̄) ./ Δt
    
    function D₁ϕ(d, q̄, q, a, Δt)
        d .= 0
        for i in eachindex(q)
            d[i,i] = - 1 / Δt
        end
    end

    function D₂ϕ(d, q̄, q, a, Δt)
        d .= 0
        for i in eachindex(q)
            d[i,i] = + 1 / Δt
        end
    end

    function Dₐϕ(d, q̄, q, a, Δt)
        d .= 0
    end


    sol = integrate(lode, HPImidpoint(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[]))
    @test relative_maximum_error(sol.q, pref.q) < 4E-4

    ref = integrate(lode, PMVImidpoint())
    @test relative_maximum_error(sol.q, ref.q) < 8*eps()


    sol = integrate(lode, HPItrapezoidal(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[]))
    @test relative_maximum_error(sol.q, pref.q) < 4E-4

    ref = integrate(lode, PMVItrapezoidal())
    @test relative_maximum_error(sol.q, ref.q) < 8*eps()

end
