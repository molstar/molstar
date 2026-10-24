# VMD: scope

What the VMD dialect does ([README](README.md)), written from the command tables by
`src/_test/scope-spec.test.ts` (`UPDATE_SPEC=1` writes it). `molid` is a molecule id or `top`; a selection applies
to the command’s molecule. Tcl: `set`, `unset`, `$name`, `[atomselect …]`, `[molinfo top]`, `puts`; atomselect
handles’ `num`, `text`, `molid`, `get`, `delete`.

## Commands

| Command | Usage | What it does in Mol\* |
| --- | --- | --- |
| `mol new` | `mol new file [type t] [waitfor all]` | Load a file given with the script, a URL, or a PDB id (`type webpdb`), with VMD’s Lines look. |
| `mol pdbload` | `mol pdbload id` | Load a PDB entry, with VMD’s Lines look. |
| `mol delete` | `mol delete molid\|all` | Remove a molecule (`all`: every one). |
| `mol on` | `mol on molid` | Show a molecule. |
| `mol off` | `mol off molid` | Hide a molecule. |
| `mol top` | `mol top molid` | Make a molecule the top one (what `top` means). |
| `mol representation` (`mol rep`) | `mol representation Style [params…]` | Set the style the next `mol addrep` takes. |
| `mol color` | `mol color Method [args]` | Set the coloring the next `mol addrep` takes. |
| `mol selection` (`mol sel`) | `mol selection {selection}` | Set the selection the next `mol addrep` takes. |
| `mol material` | `mol material Material` | Set the material the next `mol addrep` takes (Opaque, Transparent). |
| `mol addrep` | `mol addrep molid` | Add a rep to a molecule with the current style, coloring, selection and material. |
| `mol delrep` | `mol delrep repid molid` | Remove a rep of a molecule. |
| `mol modstyle` | `mol modstyle repid molid Style [params…]` | Change a rep’s style. |
| `mol modcolor` | `mol modcolor repid molid Method [args]` | Change a rep’s coloring. |
| `mol modselect` | `mol modselect repid molid {selection}` | Change a rep’s selection. |
| `mol modmaterial` | `mol modmaterial repid molid Material` | Change a rep’s material. |
| `mol showrep` | `mol showrep molid repid on\|off` | Show or hide a rep. |
| `mol list` | `mol list` | List the molecules. |
| `measure bond` | `measure bond {i j} [molid m]` | Print the distance between two atoms (VMD indices). |
| `measure angle` | `measure angle {i j k} [molid m]` | Print the angle of three atoms (VMD indices). |
| `measure dihed` | `measure dihed {i j k l} [molid m]` | Print the dihedral angle of four atoms (VMD indices). |
| `measure center` | `measure center $sel` | Print the center of an atom selection. |
| `measure rgyr` | `measure rgyr $sel` | Print the radius of gyration of an atom selection. |
| `measure minmax` | `measure minmax $sel` | Print the bounding box of an atom selection. |
| `label add` | `label add Atoms\|Bonds\|Angles\|Dihedrals molid/index …` | Label an atom, or measure a bond, angle or dihedral (atoms as molid/index). |
| `label delete` | `label delete Atoms\|Bonds\|Angles\|Dihedrals [all]` | Remove the labels or measurements of a category. |
| `color Display` | `color Display Background color` | Set the background color (`color Display Background white`). |
| `display projection` | `display projection Orthographic\|Perspective` | Orthographic or Perspective projection. |
| `display depthcue` | `display depthcue on\|off` | Fog on or off. |
| `display resetview` | `display resetview` | Reset the camera. |
| `display update` | `display update [on\|off]` | Nothing to do: Mol* draws when the scene changes. |
| `axes location` | `axes location off` | Nothing to do with `off`: Mol* shows no axes in the scene. |
| `rotate` | `rotate x\|y\|z by angle` | Rotate the view about a screen axis, in degrees (`rotate y by 90`). |
| `translate` | `translate by x y z` | Move the view by screen units (`translate by x y z`). |
| `scale` | `scale by factor` | Zoom the view (`scale by 1.2`: closer). |
| `animate goto` | `animate goto frame\|start\|end` | Show a frame of the trajectory (0-based; start, end). |
| `animate forward` | `animate forward` | Play the trajectory. |
| `animate stop` | `animate stop` | Stop playing the trajectory. |
| `render` | `render snapshot\|TachyonInternal\|… file` | Save an image (PNG; the renderer is Mol*’s). |
| `puts` | `puts text` | Print text. |
| `undo` | `undo` | Undo the last step (alone on its line; Mol*’s, VMD has none). |

Also `help [command]`, and `source file.tcl` or `play file.vmd`: a script given with this one runs in its place.

## Styles and coloring methods

Styles: `Lines`, `Bonds`, `Licorice`, `CPK`, `VDW`, `Points`, `NewCartoon`, `Cartoon`, `NewRibbons`, `Ribbons`, `Tube`, `Trace`, `QuickSurf`, `MSMS`, `Surf`, `NewSurf`. Coloring methods: `ColorID`, `Name`, `Element`, `Type`, `ResName`, `ResID`, `Chain`, `SegName`, `Molecule`, `Structure`, `Beta`, `Occupancy`, `Charge`.
Materials: `Opaque`, `Transparent` (others draw as Opaque, with a warning).

## Not supported

| Command | Why | Instead |
| --- | --- | --- |
| `proc` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `if` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `for` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `foreach` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `while` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `expr` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `switch` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `eval` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `return` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `break` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `continue` | Tcl’s control flow is not supported. | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `package` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `namespace` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `lindex` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `llength` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `lappend` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `list` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `incr` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `format` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `string` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `array` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `dict` | Tcl’s commands are not supported (set, unset, puts and $variables are). | The console runs VMD commands one by one; Tcl’s control flow and procedures are not run. |
| `measure fit` | The matrix is for $sel move, which changes coordinates. | Superpose with pymol align mobile, target. |
| `measure rmsd` | RMSD between selections is not supported yet. | pymol align prints the RMSD of a superposition. |
| `measure contacts` | Contacts are not supported. | Interactions are shown from the Structure panel. |
| `measure hbonds` | Hydrogen bonds are not supported here. | Interactions are shown from the Structure panel. |
| `measure sasa` | Surface areas are not supported. |  |
| `mol addfile` | Adding coordinates to a molecule is not supported yet. | Load a trajectory from the Load menu. |
| `mol rename` | Renaming molecules is not supported yet. |  |
| `mol fix` | Fixing molecules is not supported. |  |
| `mol free` | Fixing molecules is not supported. |  |
| `mol active` | Active and inactive molecules are not supported. | mol on and mol off show and hide one. |
| `mol inactive` | Active and inactive molecules are not supported. | mol on and mol off show and hide one. |
| `color change` | Changing VMD’s color definitions is not supported. |  |
| `color scale` | Color scales are not supported. |  |
| `color Name` | Changing the colors of names is not supported. |  |
| `color Element` | Changing the colors of elements is not supported. |  |
| `pbc` | Periodic boundaries are not supported here. | Show a simulation box from the trajectory panel. |
| `graphics` | Drawing graphics is not supported. |  |
| `draw` | Drawing graphics is not supported. |  |
| `material` | Changing materials is not supported. | mol material Opaque or Transparent. |
| `menu` | VMD’s menus are not here. |  |
| `plugin` | VMD’s plugins are not here. |  |
| `light` | Lights are set in the canvas settings. |  |
| `logfile` | Log files are not written. |  |
| `topo` | Changing the topology is not supported. | Mol* shows structures as they are loaded; it does not change coordinates. |
| `quit` | There is nothing to quit. |  |
| `exit` | There is nothing to quit. |  |
