
using GeometricIntegrators
using LinearAlgebra
using Aqua

# Aqua.test_all(GeometricIntegrators)

# ambiguities are checked for Base, LinearAlgebra, etc., which fails
# Aqua.test_all(GeometricIntegrators; ambiguities=false)

# undefined_exports doesn't seem to work well with reexport
# Aqua.test_all(GeometricIntegrators; ambiguities=false, undefined_exports=false)

Aqua.test_all(GeometricIntegrators; ambiguities = false, undefined_exports = false)
