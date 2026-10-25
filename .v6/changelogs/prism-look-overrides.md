## Prism: the look's scene defaults and overrides

- Prism: the Structure panel's sliders (was a gear) open the scene defaults in a popover: hydrogens, quality, ignore
  light, material, interior, interactions. Changing them works again (the form's values were refused by the command).
- Prism: the Animation card's Wiggle sets the wiggle's type (position, group), speed and frequency, and the tumble of
  whole copies (amplitude, speed, frequency); `animation.wiggle` takes them, each kept when absent (mode and amplitude
  too). Off stops the wiggle; the tumble stops at 0; Clear stops both.
- Prism: the Clipping section's Settings: Cut by (pixel, or whole copies by their center) for every cut, and Solid
  surface; `clipping.settings`. Cuts added later cut the same way.
- Prism: a structure's, component's or representation's inspector has its own sliders: a look over the scene defaults
  for the representations under it. Values of its own are overrides, kept when the default changes; marked by a dot on
  the sliders, by the sliders glyph on the tree rows (faint on the rows above), and listed in both popovers with Use
  default and Apply to all.
- Prism: `look.set` (the scene defaults, or a node's representations), `look.reset`, `look.get`, and their commands;
  `component.options` is gone (`look.set` runs under that id until the release).
- [Breaking] `StructureComponentManager.setOptions` moves only the representations that had the previous value of a changed
  option; one given its own value keeps it. `reset` writes options to every representation; `reprOptionOverrides`,
  `reprHasOption`, `applyReprOption`, `reprOptionValue` read and write a representation's options. Representation
  presets give their representations the material, clip objects, interior and animation of the options, as
  `addRepresentation` does.
- `Popover` in `@molstar/prism-ui`: a panel of controls from a button, which a select inside it does not close.
- Prism: a structure representation's colors offer Element (carbons by chain, as before) and Element, grey C (every
  carbon the element's grey); their swatches show their carbons, then N, O and S (were the first entries of the
  periodic table). The tree row says "grey C".
- Prism: a Shapes group in the Structure panel: each mesh loaded on its own (OBJ, PLY, VTP) by its file's name, its
  representations under it; an eye, Remove, a look (sliders) and its representation's form (color and opacity first).
  The `shape` area: `list`, `visibility`, `update`, `remove`, `focus`, and their commands. Stories: Structure panel ›
  Shapes, Shape representation (meshes made in code).
- Prism: `look.set` by node reaches volumes', particles' and shapes' representations too, and their inspectors have the
  sliders (the values they take: no hydrogens); only structure representations follow the scene defaults.
- Prism: a structure loaded is framed below the sequence band from the start (Mol*'s first framing counts in the band
  that is about to come): Reset view no longer moves it, nothing jumps when the band comes; closing or showing the band
  later glides the scene to its new frame while the camera was not moved.
