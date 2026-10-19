## Extensions: kinemage

- Prism: a Kinemage panel, a tab on the rail beside Structure (its own glyph) while a kinemage is loaded, opened when
  one loads and by "Show the Kinemage panel". Drawn as the Structure panel: Next frame in its header, then a tree of the scene tree's rows: each kinemage, its views (a click applies one), its animate frames (the row says the frame shown), its masters, and its other groups with their subgroups, each with an eye. Every control is a `kinemage` operation, so a click is a step in
  the history and Undo takes it back
- Prism: an extension's tab can be on the rail only while something holds (`visible`); `SceneRow` can be a tree of its
  own outside the scene tree (`children`, `expanded`, `depth`, `icon`) and act on a click (`onActivate`); a `kinemage`
  molecular glyph
- Storybook: Prism › Viewer › With extensions › Kinemage, a crambin kinemage (`data/examples/1crn.kin`) with its
  groups, masters, views and animate frames driven through the `kinemage` API, its commands and the Kinemage card
