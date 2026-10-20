## Graphics: interactions as cylinder impostors

- The interactions representation draws its intra-unit, inter-unit and bridge visuals as cylinder impostors, one per
  dash, like bonds; the dashed mesh stays as the fallback (`tryUseImpostor: false`, no `fragDepth`). Same groups,
  picking and marking as the mesh; geometry export now meshes the impostors as it does for bonds. Buffers: 186k atoms
  306 → 83 MB; 1M-atom assembly 0.44–0.98 GB → 0.27 GB; 5M-atom assembly (864 copies) 2.0–4.4 GB → 1.2 GB, build
  4.2 s (`auto`, quiet machine) / 20 s (`low`, loaded machine) → 3.0 s
- `physical` and `uncertainty` size themes: interaction locations get the smaller size of their two ends, as the mesh
  baked it (was the default size); a bridge gets the smallest of its four feature atoms for all its legs
- Interactions `sizeAspectRatio` (from the cylinder params) scales the mesh radius too; `colorMode: 'interpolate'`
  draws as 'default' with impostors (the mesh path ignored it)
- Fixes (perf audit 09/B1, 09/B4 and two found on the way):
  - a salt bridge hides the H-bonds between the same atoms, not itself (as documented, NGL and PLIP; ionic is off by
    default)
  - the overlap refiners also filter inter-unit contacts (the inter-unit check never matched): H-bonds within a salt
    bridge, ionic contacts within a metal coordination, hydrophobic and cation-pi contacts within pi-stacking
  - `closestHydrogenIndex` considers only hydrogens of the donor's unit (a hydrogen bonded across units gave an index
    read as an atom of the donor's unit)
  - `includeParent`: intra-unit contacts are shown or hidden by their own atoms
  - the bridge visual no longer rebuilds its geometry when an update brings a new but equal size theme
  - H-bond dashes read the donor flag per feature (it was read by member offset; latent, H-bond features come first)
