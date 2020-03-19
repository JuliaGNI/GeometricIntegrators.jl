using Weave, GeometricIntegrators

tutorial_path = "tutorial/tutorial.jmd"

# HTML
weave(tutorial_path,
  out_path="tutorial/html",
  doctype = "md2html")
