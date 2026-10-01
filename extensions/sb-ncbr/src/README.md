# SB NCBR extensions

## Partial charges

This extension visualizes partial atomic charges for atoms and residues. The extension reads charge data of a structure from a mmcif file and displays them as a color gradient on the atoms/residues. The coloring uses two gradients: one for positive charges (white-to-blue) and one for negative charges (red-to-white). The color is interpolated between the appropriate gradient based on the charge value. The extension also displays the charge values in the description label when an atom/residue is selected.

### How to use

To visualize partial charges, you need to provide a mmcif file with the structure and its charges. The charges are stored under the following categories:

```
_sb_ncbr_partial_atomic_charges_meta.id         # id of the charges (e.g. 1)
_sb_ncbr_partial_atomic_charges_meta.type       # type of the charges (optional, e.g. 'empirical')
_sb_ncbr_partial_atomic_charges_meta.method     # calculation method name (e.g. 'QEq', 'SQE+qp/Schindler 2021 (PUB_pept)')

_sb_ncbr_partial_atomic_charges.type_id         # id of the charges (pointer to _sb_ncbr_partial_atomic_charges_meta.id)
_sb_ncbr_partial_atomic_charges.atom_id         # atom id (pointer to _atom_site.id)
_sb_ncbr_partial_atomic_charges.charge          # partial atomic charge 
```
> Note that the mmcif item `_partial_atomic_charges_meta.method` is used as a description of the charge set in the UI (described in *Controls*).

The extension will automatically read the charges from the mmcif file and color the structure accordingly.

### Controls

The extension provides controls for setting the color gradient range and for selecting charge type (atom charges or residue charges).
These controls are available in Color Theme settings for 3D Representation cells in the State Tree UI.

There is also a dropdown menu for switching between charge sets.
These controls are available in Custom Model Properties settings for Model cell in the State Tree UI.
