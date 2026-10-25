# PyMOL: scope

What the PyMOL dialect does ([README](README.md)), written from the command tables by
`src/_test/scope-spec.test.ts` (`UPDATE_SPEC=1` writes it). An omitted selection is `all`, as in PyMOL.

## Commands

| Command | Usage | What it does in Mol\* |
| --- | --- | --- |
| `fetch` | `fetch code [, name [, type=cif]]` | Load an entry from the PDB (1cbs, pdb_00001cbs) or AlphaFold DB (AF-…). |
| `load` | `load filename [, object [, format=…]]` | Load a file given with the script, or a URL. |
| `delete` | `delete name` | Remove an object (a structure), or forget a named selection; `all` removes everything. |
| `enable` | `enable name` | Show an object. |
| `disable` | `disable name` | Hide an object. |
| `select` | `select [name,] selection [, enable=0]` | Name a selection and select it (one argument: the name is `sele`). |
| `deselect` | `deselect` | Clear the selection. |
| `show` | `show [representation [, selection]]` | Show atoms as a representation (no arguments: lines). |
| `hide` | `hide [representation [, selection]]` | Hide a representation of atoms (no arguments: everything). |
| `as` | `as representation [, selection]` | Show atoms as only this representation. |
| `cartoon` | `cartoon type [, selection]` | Set how a cartoon is drawn: automatic, putty or tube. |
| `color` | `color color [, selection]` | Color atoms: a PyMOL color name, 0xRRGGBB or [r, g, b]; `atomic` colors by element. |
| `util.cbag` | `util.cbag [selection]` | Color carbons green and the other atoms by element. |
| `util.cbac` | `util.cbac [selection]` | Color carbons cyan and the other atoms by element. |
| `util.cbam` | `util.cbam [selection]` | Color carbons magenta and the other atoms by element. |
| `util.cbay` | `util.cbay [selection]` | Color carbons yellow and the other atoms by element. |
| `util.cbas` | `util.cbas [selection]` | Color carbons salmon and the other atoms by element. |
| `util.cbaw` | `util.cbaw [selection]` | Color carbons white and the other atoms by element. |
| `util.cbab` | `util.cbab [selection]` | Color carbons slate and the other atoms by element. |
| `util.cbao` | `util.cbao [selection]` | Color carbons orange and the other atoms by element. |
| `util.cbap` | `util.cbap [selection]` | Color carbons purple and the other atoms by element. |
| `util.cbak` | `util.cbak [selection]` | Color carbons pink and the other atoms by element. |
| `util.cnc` | `util.cnc [selection]` | Color the atoms other than carbon by element. |
| `util.cbc` | `util.cbc [selection]` | Color by chain. |
| `util.chainbow` | `util.chainbow [selection]` | Color each chain as a rainbow. |
| `util.rainbow` | `util.rainbow [selection]` | Color as a rainbow from the N to the C terminus. |
| `util.ss` | `util.ss [selection]` | Color by secondary structure (Mol*’s colors). |
| `spectrum` | `spectrum [expression [, palette [, selection]]]` | Color along a property: count or resi (residue number), b (B-factor), chain. |
| `bg_color` | `bg_color color` | Set the background color. |
| `set` | `set name, value [, selection]` | Change a setting (the ones Mol* has: see `help set`). |
| `unset` | `unset name [, selection]` | Put a setting back to its default. |
| `zoom` | `zoom [selection [, buffer [, state [, complete [, animate]]]]]` | Frame a selection (no selection: everything). |
| `orient` | `orient [selection [, state [, animate]]]` | Orient along the principal axes (with a selection: its structure’s, then framed on it). |
| `center` | `center [selection [, state [, origin [, animate]]]]` | Center the view on a selection, keeping the zoom. |
| `reset` | `reset` | Reset the camera. |
| `turn` | `turn axis, angle` | Rotate the camera about an axis of the screen, in degrees. |
| `move` | `move axis, distance` | Move the camera along an axis of the screen, in Ångströms. |
| `distance` (`dist`) | `distance [name,] selection1, selection2 [, cutoff [, mode]]` | Measure the distance between two atoms, or every pair of two sets within a cutoff. |
| `angle` | `angle [name,] selection1, selection2, selection3` | Measure the angle of three atoms. |
| `dihedral` | `dihedral [name,] selection1, selection2, selection3, selection4` | Measure the dihedral angle of four atoms. |
| `label` | `label [selection [, expression]]` | Label atoms: resn, resi, chain, name, elem, b, quoted text, + and "%s" % (…); empty clears. |
| `align` | `align mobile, target` | Superpose by sequence alignment (Mol*’s chain superposition); prints the RMSD. |
| `super` | `super mobile, target` | Superpose by structure (TM-align); prints the RMSD and TM-score. |
| `cealign` | `cealign target, mobile` | Superpose by structure (TM-align, Mol*’s structure alignment); prints the RMSD. |
| `pair_fit` | `pair_fit mobile1, target1 [, mobile2, target2, …]` | Superpose by pairs of atoms (mobile, target, mobile, target, …). |
| `remove` | `remove selection` | Leave whole entities or chains out of the scene (hydrogens: `remove hydro`). |
| `png` | `png filename [, width [, height [, dpi [, ray]]]]` | Save an image (PNG). |
| `ray` | `ray [width [, height]]` | Nothing to do: Mol* renders live (`png` saves the image). |
| `draw` | `draw [width [, height]]` | Nothing to do: Mol* renders live (`png` saves the image). |
| `count_atoms` | `count_atoms [selection]` | Count the atoms of a selection. |
| `get_names` | `get_names` | List the objects and the named selections. |
| `get_object_list` | `get_object_list` | List the objects. |
| `get_view` | `get_view` | Describe the camera. |
| `frame` | `frame number` | Show a frame of the trajectory (1-based, as PyMOL’s states). |
| `mplay` | `mplay` | Play the trajectory. |
| `mstop` | `mstop` | Stop playing the trajectory. |
| `undo` | `undo` | Undo the last step (alone on its line). |
| `redo` | `redo` | Redo what Undo took back. |

Also `help [command]`, and `@file.pml` or `load file.pml`: a script given with this one runs in its place.

## Representations

| PyMOL | Mol\* |
| --- | --- |
| `cartoon` | `cartoon` |
| `ribbon` | `backbone` |
| `sticks` | `ball-and-stick` (sizeAspectRatio: 1) |
| `lines` | `line` |
| `spheres` | `spacefill` |
| `surface` | `molecular-surface` |
| `mesh` | `molecular-surface` (visuals: ["molecular-surface-wireframe"]) |
| `dots` | `point` |
| `ellipsoids` | `ellipsoid` |
| `nb_spheres` | `spacefill`, atoms without bonds |
| `nonbonded` | `point`, atoms without bonds |
| `putty` | `putty` |
| `everything` | every type |
| `labels` | not supported: Labels come from `label selection, expression`. |
| `slice` | not supported: Map slices are not shown; maps are phase 2. |
| `extent` | not supported: Map extents are not shown. |
| `cgo` | not supported: CGO objects are not read. |
| `cell` | not supported: Show the unit cell from the Structure panel. |
| `callback` | not supported: Callbacks are Python. |
| `dashes` | not supported: Measurements are made with `distance`. |
| `angles` | not supported: Measurements are made with `angle`. |
| `dihedrals` | not supported: Measurements are made with `dihedral`. |

## Settings

| Setting | What | Default |
| --- | --- | --- |
| `bg_rgb` | The background color. | `black` |
| `bg_color` | The background color. | `black` |
| `orthoscopic` | Orthographic projection (1) or perspective (0). | `0` |
| `field_of_view` | The camera’s field of view, in degrees. | `20` |
| `depth_cue` | Fog: far atoms fade (1) or not (0). | `1` |
| `fog` | Fog: far atoms fade (1) or not (0). | `1` |
| `transparency` | Transparency of surfaces, 0 (opaque) to 1. | `0` |
| `cartoon_transparency` | Transparency of cartoons, 0 to 1. | `0` |
| `stick_transparency` | Transparency of sticks, 0 to 1. | `0` |
| `sphere_transparency` | Transparency of spheres, 0 to 1. | `0` |
| `sphere_scale` | The size of spheres (PyMOL’s 1 is Mol*’s). | `1` |
| `stick_radius` | The radius of sticks (PyMOL’s 0.25 is Mol*’s). | `0.25` |
| `ray_opaque_background` | Whether png keeps the background (1) or makes it transparent (0). | `1` |

Any other setting fails, naming these.

## Not supported

| Command | Why | Instead |
| --- | --- | --- |
| `python` | Python blocks are not run. | Scripts here are PyMOL commands; Python is not run. |
| `run` | Python scripts are not run. | Scripts here are PyMOL commands; Python is not run. |
| `cmd` | PyMOL’s Python API is not run. | Scripts here are PyMOL commands; Python is not run. |
| `alter` | Changing atom properties with Python expressions is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `alter_state` | Changing coordinates is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `iterate` | Reading properties with Python is not supported yet. | `count_atoms` counts a selection. |
| `iterate_state` | Reading coordinates with Python is not supported. | Scripts here are PyMOL commands; Python is not run. |
| `create` | New objects from selections are not supported. | Mol* shows parts of a structure as components: `show sticks, …`, `select name, …`. |
| `extract` | New objects from selections are not supported. | Mol* shows parts of a structure as components: `show sticks, …`, `select name, …`. |
| `copy` | Copying objects is not supported. | Load it again with `fetch` or `load`. |
| `split_states` | Splitting states into objects is not supported. | `frame n` shows one state. |
| `fab` | Building peptides is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `fuse` | Joining objects is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `h_add` | Adding hydrogens is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `h_fill` | Adding hydrogens is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `mutagenesis` | The mutagenesis wizard is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `wizard` | Wizards are not supported. |  |
| `edit` | Editing is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `bond` | Adding bonds is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `unbond` | Removing bonds is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `sort` | Sorting atoms is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `sculpt_activate` | Sculpting is not supported. | Mol* shows structures as they are loaded; it does not edit them. |
| `mset` | Movies are not supported. | Use the snapshots and animation panels. |
| `mview` | Movies are not supported. | Use the snapshots and animation panels. |
| `mdo` | Movies are not supported. | Use the snapshots and animation panels. |
| `save` | Saving files from a script is not supported. | Save the session from the menu (⌘S). |
| `symexp` | Symmetry mates as objects are not supported. | Show symmetry from the Structure panel. |
| `ramp_new` | Color ramps are not supported. | `spectrum` colors along a property. |
| `set_name` | Renaming objects is not supported yet. |  |
| `order` | Ordering objects is not supported. |  |
| `group` | Grouping objects is not supported. |  |
| `ungroup` | Grouping objects is not supported. |  |
| `dss` | Mol* assigns secondary structure when it loads; there is nothing to recompute. |  |
| `set_view` | Setting PyMOL’s view matrix is not supported yet. | `zoom`, `orient` and `turn` move the camera. |
| `scene` | Scenes are not supported yet. | Use the snapshots panel. |
| `isomesh` | Maps are not supported yet. | Load density from the Volume panel. |
| `isosurface` | Maps are not supported yet. | Load density from the Volume panel. |
| `map_new` | Maps are not supported. |  |
| `quit` | There is nothing to quit. |  |
