using GeometricIntegrators
using ExplicitImports
using Test

@testset "ExplicitImports" begin
    test_explicit_imports(GeometricIntegrators;
        # the package brings its dependencies in with `using`, and re-exports them
        no_implicit_imports = false,
        # several imported names are internal to GeometricIntegratorsBase and its siblings
        all_explicit_imports_are_public = false,
        # the source reaches internal names of dependencies by qualified access
        all_qualified_accesses_are_public = false)
end
