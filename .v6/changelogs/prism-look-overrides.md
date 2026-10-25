## Prism: the look's scene defaults and overrides

- Prism: the Structure panel's sliders (was a gear) open the scene defaults in a popover: hydrogens, quality, ignore
  light, material, solid surface, interior (the clip objects, interactions and animation there until they move to their
  tabs). Changing them works again (the form's values were refused by the command).
- Prism: a structure's, component's or representation's inspector has its own sliders: a look over the scene defaults
  for the representations under it. Values of its own are overrides, kept when the default changes; marked by a dot on
  the sliders, by the sliders glyph on the tree rows (faint on the rows above), and listed in both popovers with Use
  default and Apply to all.
- Prism: `look.set` (the scene defaults, or a node's representations), `look.reset`, `look.get`, and their commands;
  `component.options` is gone (`look.set` runs under that id until the release).
- `StructureComponentManager.setOptions` moves only the representations that had the previous value of a changed
  option; one given its own value keeps it. `reset` writes options to every representation; `reprOptionOverrides`,
  `reprHasOption`, `applyReprOption`, `reprOptionValue` read and write a representation's options. Representation
  presets give their representations the material, clip objects, interior and animation of the options, as
  `addRepresentation` does.
- `Popover` in `@molstar/prism-ui`: a panel of controls from a button, which a select inside it does not close.
