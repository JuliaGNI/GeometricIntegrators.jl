using GeometricIntegrators
using Weave

ENV["GKSwstype"] = "100"

tutorial_path = joinpath("tutorial", "tutorial.jmd")
build_path = joinpath("src", "tutorial")

weave(tutorial_path,
  out_path=build_path,
  doctype = "github")
