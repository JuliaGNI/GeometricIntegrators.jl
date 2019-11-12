using Weave, GeometricIntegrators

tutorial_path = "src/tutorial.jmd"

# HTML
weave(tutorial_path,
  out_path="weave/html",
  doctype = "md2html")
