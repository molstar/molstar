# Command language: to do

What is missing, and what was considered and deferred; an entry leaves when it is done (its `_spec` section says so
then). The spec: [\_spec](_spec/README.md).

## Next

- **Readable targets**: an atom set written back as an address list (`A`, `B:12-14`) instead of source-index ranges
  (`compactExpression`), so the inspector reads it and MVS takes it as schema selectors.
- **A key for the console** (\` was proposed): Prism's keybindings name commands; `script.console` could take one.

## Later

- **MolViewSpec export** of a script's result: in `packages/mvs`, separately.
- **A tool for the agent** (`run_script` in a skill), for "run this .pml".
- **Per-representation colors and sizes in PyMOL** (`cartoon_color`, `stick_radius` on one type only): the style layers
  are per atom over every type.
- **PyMOL**: `set_view`, `scene`, maps (`isomesh`), `iterate` with `print`, `set_name`. **VMD**: `mol addfile`,
  `atomselect macro`, `mol rename`, `HBonds`.
- **Jmol or ChimeraX** commands: another table on the same IR.
