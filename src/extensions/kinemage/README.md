# Kinemage extension

This extension adds support for the Kinemage molecular graphics format based on the
[kinemage](http://kinemage.biochem.duke.edu/static/files/PDFs/format-kinemage.pdf) format specification.

It currently supports the following features:
- Drag-and-drop of Kinemage files into the display area
- Display of @ball, @sphere, @vector, @dot, @ribbon, and @triangle objects
- Hovering over objects to see their labels (if present)

Currently unsupported features include:
- master and submaster selections of visible objects
- Setting one or more viewpoints
- animations
- @label and @ring objects
- @hsvcolor keyword for coloring by hue, saturation, and value
- 'fore' and 'rear' keywords for different front and back colors

Current limitations include:
- Lines and triangles are a single color, not colored by vertex
- Line segments in Mol* do not support end-caps for wide lines, so there are artifacts in highly-curved lines
- The default perspective view and white background for Mol* differs from that of Kinemage