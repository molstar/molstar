## Particles

- [Breaking] Particles are opt-in: `ParticlesRegistry` (`@molstar/plugin/particles-registry`) holds the particle
  formats, transformers, representations, themes and animation, and `DefaultRegistry` none of them; the classic viewer
  no longer reads particle files.
- Prism: a lazy particles feature, loaded by a particle file (`.star`, `.tbl`, `.ndjson`, `.em`, `.simularium`), an
  mmCIF holding particles (CellPack, PetWorld; by its content) or a session with particles; lazy features can claim
  files by content (`owns.sniff`).
- Small real particle files for tests and stories in `data/examples/particles/` (CryoET Data Portal, RELION, Dynamo,
  Artiatomi, Simularium), with their sources and licenses.
