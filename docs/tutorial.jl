using Weave, GeometricIntegrators

tutorial_path = joinpath(dirname(pathof(GeometricIntegrators)), "../docs/src", "tutorial.jmd")

# Markdown
weave(tutorial_path,
  out_path="src/weave",
  doctype="pandoc")
