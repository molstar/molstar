## Particles: particle packs

- A particle pack (`.mspk`, `formats/particles/pack.ts`): a particle list and its targets' sites as raw arrays read as
  views of the file, the targets as ranges. The yeast cell (31.5 M particles) in Node: the list in 0.1 s instead of
  10.3 s, 1.96 GB instead of 6.99 GB.
- Loaded in the Viewer and Prism (`particle_pack`, `.mspk`): the yeast cell's list in 0.17 s instead of 9.5 s, its heap
  3.77 GB instead of 7.29 GB, the first frames 7.7 s sooner.
- `particle-pack` (`cli/particle-pack`): mmCIF particle assemblies (cellPACK, PetWorld; `.cif`, `.bcif`) to packs, and
  packs back to mmCIF.
- [Breaking] `ParticleList.targets`, `entities`, `compartments` and `radii` are arrays or runs (`ParticleValues`), `keys`
  optional: a pack's list holds a target's values once (the yeast cell's heap in the Viewer 3.77 → 3.17 GB).
- Particle packs version 2: a target's sites quantized (8 bytes a site, within 0.01 Å; `ParticleSites.quantized`, read
  a site at a time by `ParticleSites.getPosition` and `getRadius`): the yeast cell's pack 1.25 → 1.09 GB.
- A particle pack's sites drawn as they are (`Spheres.fromSites`: quantized centers, 8 bytes a site on the GPU, no float
  copy; WebGL2 and WebGPU): the yeast cell's scene 1.80 -> 1.50 GB of JS heap, its GPU data 531 -> 376 MB, its
  representations 6.6 -> 5.9 s (WebGPU).
- A pack's sites sized by their types (`siteType` sizes: a table of the radii by type, not a size a site) where a size
  theme sizes them by their radii: the yeast cell's GPU data 453 -> 376 MB and its CPU layers 228 -> 150 MB (WebGPU).
- Particle packs version 6: tile-ordered packs (`particle-pack --tiles 16`, `writeParticlePack(list, { tiles })`; a pack
  written again with `particle-pack in.mspk out.mspk`): the particles by tile in Morton order and a tile directory, a
  tile set of that size or coarser made without a pass over the particles (the yeast cell at 16³ in Chrome: 3.8 s ->
  0.4 s); a target's reach and volume stored, a tile set's catalog without a pass over the sites; the tiles'
  opacities stored, a view by the camera leaving out the tiles hidden without computing them (the yeast's fill at 16³
  9.1 -> 4.0 s). A view's `opticalDepth`: the occlusion's threshold, lower leaving out more (faster, at times a gap).
- A tile-ordered particle pack drawn whole (not as tiles) is drawn by its tiles on WebGPU: the yeast cell in 3.7 s
  instead of 5.3 s (a pack written by target 2.2 s), so one file serves both.
