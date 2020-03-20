using GeometricIntegrators
using Weave

tutorial_path = joinpath(dirname(pathof(GeometricIntegrators)), "../docs/tutorial", "tutorial.jmd")
build_path = joinpath(dirname(pathof(GeometricIntegrators)), "../docs/src/tutorial")

weave(tutorial_path,
  out_path=build_path,
  doctype = "github")
