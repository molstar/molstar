# Kinemage extension

This extension adds support for the Kinemage molecular graphics format based on the
[kinemage format specification](http://kinemage.biochem.duke.edu/static/files/PDFs/format-kinemage.pdf).

It currently supports the following features:
- Drag-and-drop of Kinemage files into the display area
- Open File can open Kinemage files from the local filesystem
- Display of @ball, @sphere, @vector, @dot, @ribbon, and @triangle lists
- Coloring of objects by vertex color, or by a single color for the entire list
- Hovering over objects to see their labels (if present)
- When there are views defined, controls are added to the right panel; when selected, they shift the view
- Control panel names are based on the @pdbfile or @caption in the Kinemage file if there is one
- Lines are split in half, with each half colored by and labeled by the nearest vertex
- Master and submaster selections of visible objects
- Group and subgroup hierarchy with buttons to control visibility
- animate/2animate: First entry turned on to start, changing visibility of Animate button cycles through them

Currently unsupported features include:
- @pointmaster lists controlling visibility
- @label and @ring lists
- @hsvcolor keyword for coloring by hue, saturation, and value
- 'fore' and 'rear' keywords for different front and back colors

Current limitations include:
- Triangles are a single color, not colored by vertex (Mol* does not support per-vertex coloring for these primitives)
- Line segments in Mol* do not support end-caps for wide lines, so there are artifacts in highly-curved lines
- The default perspective view and white background for Mol* differs from that of Kinemage (though selecting a view from the
  State Tree will switch it to orthographic)