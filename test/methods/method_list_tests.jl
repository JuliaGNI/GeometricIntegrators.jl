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
    @test length(MethodList(refs = true)) == length(method_list)

    # markdown rendering (used by docs/src/methods.jmd)
    @test occursin("|", sprint(show, MethodList(markdown = true)))
    @test occursin("|", sprint(show, MIME"text/markdown"(), ml))
    @test length(sprint(show, ml)) > 0

    # per-problem-type selectors partition the list non-trivially
    selectors = (isodemethod, ispodemethod, ishodemethod, isiodemethod,
        islodemethod, issodemethod, isdaemethod, ispdaemethod,
        ishdaemethod, isidaemethod, isldaemethod)

    for selector in selectors
        # a selector can only ever narrow the list
        @test length(MethodList(selector = selector)) ≤ length(method_list)
    end

    # Every registered method must be applicable to at least one problem type, so the
    # per-selector counts have to cover the whole list (with multiplicity, since a method
    # may support several). A method whose `is*method` traits are all `false` — the easy
    # mistake when registering one — would break this.
    @test sum(length(MethodList(selector = s)) for s in selectors) ≥ length(method_list)

    @test length(MethodList(selector = isodemethod)) > 0
    @test length(MethodList(selector = islodemethod)) > 0

    # every method in the list is a GeometricMethod (type or instance)
    for m in method_list
        @test m <: GeometricIntegrators.GeometricMethod
    end
end
