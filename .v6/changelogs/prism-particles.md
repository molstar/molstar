## Particles

- [Breaking] Particles are opt-in: `ParticlesRegistry` (`@molstar/plugin/particles-registry`) holds the particle
  formats, transformers, representations, themes and animation, and `DefaultRegistry` none of them; the classic viewer
  composes it and reads particle files as before.
- Prism: a kept scene with particles comes back after a reload; a downloaded or gzipped mmCIF holding particles is
  read as particles; a particle's plain-text label has no markup.
- Prism: a lazy particles feature, loaded by a particle file (`.star`, `.tbl`, `.ndjson`, `.em`, `.simularium`), an
  mmCIF holding particles (CellPack, PetWorld; by its content) or a session with particles; lazy features can claim
  files by content (`owns.sniff`).
- Small real particle files for tests and stories in `data/examples/particles/` (CryoET Data Portal, RELION, Dynamo,
  Artiatomi, Simularium), with their sources and licenses.
- Particles without radii are drawn at half the typical distance between neighbors (the `particle-size` theme's
  `defaultSize` is 0, "from the spacing", instead of 1 Å; `getParticleSpacingRadius`), so cryo-ET picks show at about
  their size.
- Prism: the canvas' trajectory bar follows a particle trajectory's frames (Simularium) as a structure's, through
  `particles.frame`, `.play` and `.stop`; a particle list's frame row steps it back as a structure's does.
- A particle reads "Particle 165 · tomo 1 · class 2 · cc 0.4213" in bold, its list's name beneath (was "list | #165 |
  tomo 1 | …"); several are "N particles".
- Prism: a `.molx` session holding particles loads the particles feature first; a list (or map) loaded from an address
  is named by its file, not the whole URL; long names wrap in the hover label and row tooltips instead of overflowing.
