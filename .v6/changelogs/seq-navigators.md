## Seq* (next): navigators

- A `navigator` spec node bound to sections, in four forms with builders (`stripNavigator`, `miniatureNavigator`,
  `ideogramNavigator`, `overviewNavigator`), and a `bands` layer type (Giemsa stains, centromere triangles).
- The ruler rail on every linear section by default (`navigator.form: rail`), a slim rail without a ruler, and a rows
  rail beside an alignment's rows: drag to pan, grips to zoom, click to centre, sweep to zoom to a region,
  double-click to fit; one gesture, one history step.
- `mountNavigator` places a navigator part beside the view; an export includes navigators on request
  (`export.navigators` on the view, or `renderSvg`'s `navigators` option).
