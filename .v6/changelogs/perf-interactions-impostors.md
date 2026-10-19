## Graphics: interactions as cylinder impostors

- The interactions representation draws its intra-unit, inter-unit and bridge visuals as cylinder impostors, one per
  dash, like bonds; the dashed mesh stays as the fallback (`tryUseImpostor: false`, no `fragDepth`, export). Same
  groups, picking and marking as the mesh. Buffers: 186k atoms 306 → 83 MB; 1M-atom assembly 0.44–0.98 GB → 0.27 GB;
  5M-atom assembly (864 copies) 2.0–4.4 GB → 1.2 GB, build 4.2 s (quiet machine) / 20 s (`low`) → 3.0 s
- `physical` size theme: interaction and bridge locations get the smaller vdW radius of their ends (was the default
  size), so impostors are as thick as the mesh bakes them
- Interactions `colorMode: 'interpolate'` draws as 'default' with impostors (the mesh path ignored it)
