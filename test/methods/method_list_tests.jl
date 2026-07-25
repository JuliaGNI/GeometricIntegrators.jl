using GeometricIntegrators
using GeometricIntegrators.Integrators: MethodList, method_list
using Test


# `MethodList` builds the documentation tables in `docs/src/methods.jmd` by calling
# `order(m)`, `isexplicit(m)`, `issymmetric(m)`, `issymplectic(m)` and the eleven
# `is*method(m)` predicates on every registered method *type*. Registering a new
# method without the corresponding trait overloads therefore breaks the docs build
# rather than the test suite — which is exactly what this testset guards against.
#
# (This replaces the former `test/methods/test_list.jl`, which was never run and
# called `Methods.MethodList()`; the `Methods` submodule was merged into
# `Integrators`, and `MethodList` is not exported.)

@testset "$(rpad("Method list",80))" begin

    # every registered method must survive the full trait sweep
    ml = MethodList()
    @test length(ml) == length(method_list)
    @test length(ml) > 0

    # the `refs = true` variant is the one the docs use; it must not throw
    @test length(MethodList(refs=true)) == length(method_list)

    # markdown rendering (used by docs/src/methods.jmd)
    @test occursin("|", sprint(show, MethodList(markdown=true)))
    @test occursin("|", sprint(show, MIME"text/markdown"(), ml))
    @test length(sprint(show, ml)) > 0

    # per-problem-type selectors partition the list non-trivially
    for selector in (isodemethod, ispodemethod, ishodemethod, isiodemethod,
                     islodemethod, issodemethod, isdaemethod, ispdaemethod,
                     ishdaemethod, isidaemethod, isldaemethod)
        n = length(MethodList(selector=selector))
        @test 0 ≤ n ≤ length(method_list)
    end

    @test length(MethodList(selector=isodemethod)) > 0
    @test length(MethodList(selector=islodemethod)) > 0

    # every method in the list is a GeometricMethod (type or instance)
    for m in method_list
        @test m <: GeometricIntegrators.GeometricMethod
    end

end
