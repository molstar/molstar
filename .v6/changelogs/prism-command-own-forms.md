## Prism

- Add representation… and Set coloring… (a component's inspector, and Show as… / Color by theme… in the command menu)
  open their own forms: the representation types the component supports as icon rows (those it shows marked), with
  the picked type's size, quality and opacity; the color themes that apply, painted with their colors on the structure,
  and the uniform color's swatches. The themes are grouped by the registry's categories, as in the inspector's theme
  select; a click on a theme (or on one of the uniform color's swatches) applies it. Typing filters, the arrows move the
  pick, Enter runs.
- Fixed: Color by theme… opened on its query-language field (PyMOL, VMD, …) instead of asking for the theme.
- prism-ui: a command can have its own form in the command menu (`CommandMenuProps.forms`, `CommandFormView`), which
  sets form-only values with `CommandInvocation.setParam`.
- `representation.color` with `theme: 'uniform'` takes `color` as the theme's color.
- A component's representations in its inspector have their row's ⋯ menu (Remove, Type, Color, Edit…). Their eye and ⋯
  sit left of the chevron and show on the row's hover or focus (always on a touch screen), as the tree's do.
- Fixed: a row's ⋯ menu ran off the window when an entry opened its sub-entries; it moves to stay in view, and scrolls
  when taller than the window.
- A command's query language (PyMOL, …) shows in its sentence only while the target is typed as a query; Show as… and
  Color by theme… read as sentences ("Polymer in 1TQN as Cartoon", "ligand in 1TQN by Chain Id").
