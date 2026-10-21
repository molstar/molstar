# Prism particles: plan (working notes, deleted before the PR)

This comes before the clipping extension (`.v6/clipping-extension/plan.md`), which needs particles as a clippable kind of node. Before the PR, what holds moves into `packages/prism/viewer/_spec/`, the TODO and the docs, and this folder is deleted.

## What exists (survey, 2026-10-21)

**Model and IO** (`packages/model/src/model/particles/`, `formats/particles/`):
- `ParticleList` holds:
  - `coordinates`, plus optional `rotations` (quaternions), `radii` and `cell`;
  - `fibers`, in CSR form;
  - `attributes` (cc, class, score), with min and max;
  - entities and compartments;
  - `targets` (an id per particle), `targetInfo`, and an optional `targetMapping` from an id to a structure, shape or volume.
- `ParticleTrajectory` is for Simularium.
- Readers:
  - RELION STAR, with tomogram and micrograph filters, the pixel size and a radius;
  - Dynamo TBL;
  - CryoET Data Portal NDJSON;
  - Artiatomi `.em`;
  - mmCIF (CellPack and PetWorld assemblies);
  - Simularium.
- **Every test builds its data inline. No particle file is committed.**

**Plugin** (`packages/plugin/core`):
- State objects: `SO.Particle.{Trajectory, List, Representation3D}`.
- Format providers: `state/formats/particles/` (six of them, in the default registry).
- Presets: `simpleVisuals` (spacefill plus unit cell) and `complexVisuals` (targets, spacefill without targets, fibers, and a color theme from what exists).
- Transforms:
  - `ParticleListFromTrajectory`, with `frameIndex`;
  - `ParticleListWithTargets`, a decorator that takes files, URLs or a PDB id;
  - `ParticlesRepresentation3D`, with `{ type, colorTheme, sizeTheme }`. Its clip is at `type.params.clip`;
  - `ParticleListUnitcell3D`, a Shape representation.
- The hierarchy is `managers.particles` (`ParticleListRef.representations`). The unit cell and the trajectory are not in it.

**Graphics:**
- Representations: `spacefill`, `orientation`, `fibers`, `target` (structure, volume, shape and blob-surface sub-representations).
- Color themes: particle-index, -entity, -compartment, -hierarchy, -attribute.
- Size theme: particle-size.

**Prism** (`packages/prism/viewer`):
- **Loading:** particle files load through `load.files` (the default registry). There is no particle scope in the loaders.
- **API:** `particles.{list, representation, update, visibility, remove}` (`api/particles.ts`), registered inside `VolumeFeature` (`features/volume/index.ts:9-24`), with commands in `commands/volume.ts:53-70` and tests in `api/_test/particles.test.ts`.
- **Scene tree:**
  - The kinds `particles` and `particles-representation` exist (`model/scene.ts:25-26`, `sources/scene.ts:755-771`).
  - The inspector shows only Particles and Targets, plus the raw param form for a representation.
- **Node selectors have no particle kinds** (`api/nodes.ts:22`, `nodeKind()`). This is listed as a limit in `_spec/nodes.md:32` and `TODO.md:63`.
- **Stories:** one, `data.stories.ts` `ParticlesStar`, with a 400-pick STAR generated in code. The mock scene (`testing/scene.ts`) has no particle nodes.

## P0. Particles are opt-in: one registry entry, nothing by default (David, 2026-10-21)

**Today particles are in everything:**
- `DefaultRegistry` (`plugin/core/src/default-registry.ts`) spreads particles over several entries:
  - `DefaultThemes`, its `particles` scope;
  - `DefaultParticleRepresentations`;
  - `DefaultActions`: six format transforms, `ParticleListFromTrajectory`, `ParticleListWithTargets`, `ParticleListUnitcell3D`, `ParticlesRepresentation3D`;
  - `DefaultFormats` (`BuiltInParticlesFormats`);
  - `DefaultAnimations` (`AnimateParticleTrajectory`).
- `state/transforms/catalog.ts:16-21, 48-50` imports the particle formats and transforms for their side effects.
- Prism takes `DefaultRegistry` as is (`sources/base.ts:25`).

**Change:**
- **One entry, `ParticlesRegistry`** (`@molstar/plugin/registry/particles`), holding the formats, transforms, representations, the `particles` theme scope and the trajectory animation.
- **`DefaultRegistry` and the transform catalog no longer import any of it.** Particle transforms are registered by importing the entry, which is what a session with particles needs. A session holding particles without the entry fails, naming the missing transformer, as any unknown transformer does.
- **Apps opt in.** The classic `apps/viewer` composes `ParticlesRegistry`, so the classic viewer keeps particles. The default `PluginSpec` does not.
- **Prism:** `ParticlesFeature` is a *lazy* built-in feature (`LazyViewerExtension`, `builtIn`). Its stub declares `owns: { formats, fileExtensions: ['star', 'tbl', 'ndjson', 'em', 'simularium'], transformers, representations, themes }`, so a particle file or session loads its chunk first (`spec.ts:480`). Its `plugin.registry` is `[ParticlesRegistry]`.
  - mmCIF particles (`mmcif_particles`, sniffed from `.cif`/`.bcif`) cannot be owned by extension. It is decided in step 2: a sniff in the stub, or kept with the structure loaders.
- **Recorded in** `breaking-v6-changes.md`, which covers `DefaultRegistry` without particles and the catalog not registering them. The `default-registry` tests change.
- **Left for the `PluginContext` rewrite** (Prism step 2, see memory): `PluginContext` still builds `representation.particles` and `managers.particles.hierarchy` unconditionally, and imports their code. Moving them into the entry is part of that rewrite; noted in the TODO.

## Gaps, and what to do

| #   | Gap | Work |
| --- | --- | --- |
| P1  | No particle kinds in node selectors | `particles` and `particles-representation` in `NodeKinds` and `nodeKind()`, with attributes (`label`, `type`, `format`). Update `nodes.md` and the TODO. The clipping extension's `node` target needs this. |
| P2  | Particles live inside `VolumeFeature` | Their own lazy `features/particles/` (operations, commands, the scene parts, the stub's `owns`; see P0), and a `viewer-spec.md` section. |
| P3  | Thin inspector | List facts: count, format, targets, entities, compartments, attribute ranges, the cell or bounds. A representation's row says its type and color theme, as structure rows do. Quick settings: color by entity, compartment, hierarchy or attribute (the attribute picked), size and radius. Volume's inspector is the model. |
| P4  | No row actions | Row menu: Add representation (spacefill, orientation, fibers, target), focus or camera to the list, highlight on hover, Remove. |
| P5  | The unit cell is invisible in the tree | Show it as a child row of the list (eye, remove), as a structure's unit cell shows. |
| P6  | No targets in Prism | `particles.targets`: `ParticleListWithTargets` by a loaded file, a URL or a PDB id, per target id, with `includeParent`. Inspector: the targets table, each assigned or not. |
| P7  | Loading choices | A particle scope in the loaders (format, and RELION's tomogram or micrograph filter, which defaults to the first tomogram). The pixel size and radius are set on load and editable in the inspector. |
| P8  | Simularium frames | The trajectory row (`views/structure/trajectory-row.tsx`) for a particle trajectory: frame, scrub and play. |
| P9  | Tests and stories | Particle nodes in the mock scene. Stories per format and for the tree, the inspector, targets and frames. |
| P10 | Hierarchy rebuild cost (`TODO.md:161`) | Measure with a 100k-particle list when we are there. It is fixed only if it shows. |

Clip is not here: the clipping extension covers it, through `type.params.clip`, once P1 lands.

## Fixtures

**`data/examples/particles/`**, small real files, each with its source and license in a `README.md`:

| File | Source | Size |
| --- | --- | --- |
| `ribosomes.ndjson` | CZ CryoET Data Portal, a ribosome `orientedPoint` annotation | trimmed to about 200 picks, under 100 KB |
| `hiv-capsid.star` | RELION subtomogram tutorial, EMPIAR-10164 | trimmed to one or two tomograms, under 100 KB |
| `picks.tbl` | Dynamo tutorial table | a few hundred rows |
| `motivelist.em` | Artiatomi | a few hundred particles; generated by our own writer if no small public one exists |
| `simularium-small.simularium` | Simularium example | a few frames, under 500 KB |

**Targets:** pair a list with an existing structure (`1crn.cif`, `1cbs_updated.cif`), or a small map if one is needed. Files are read from disk by tests and served by `ExampleFiles` in stories, so nothing reaches the network.

**Generators stay** for the size sweeps (P10) and for inline quick tests.

## Steps

**Done (2026-10-21):** fixtures (step 1), P0 (step 2) and P2 (step 3):
- `ParticlesRegistry` is classified as a default-composition module in `scripts/workspace/catalogs.json`.
- The stub has `owns.sniff` (`ExtensionNeed.files` in `ensureFor`).
- The particle area is now `features/particles/operations.ts`, as the boundary rule wants.
- P1 is done: `particles` and `shape` node kinds.
- P3–P5 are done:
  - three new operations: `particles.color`, `particles.focus` and `particles.unit-cell`;
  - rows say what they show (theme, swatch, targets) and have menus (Add representation, Color by, Focus, unit cell);
  - the unit cell is a row of its own;
  - the inspector shows the list's facts (source, orientations, targets, entities, compartments, fibers, attribute ranges, cell).
- P6 is done:
  - `particles.targets`: one target per call, a file or a PDB id or URL (network: asks for a URL), adding or updating `ParticleListWithTargets`, and a `target` representation when there is none;
  - the list's Targets fact names what stands for each;
  - the row menu opens the command;
  - operations and the scene read the list at the end of its decorators (`getDecoratorRoot`).
- Not handled: a `.molx` session holding particles does not load the feature. `session.open` reads transformers only from `.molj`. To do with P9.


1. Fixtures (needs your OK to download; see the questions).
2. P0, `ParticlesRegistry` and taking particles out of `DefaultRegistry` and the catalog. The classic viewer composes it. Tests show that a plugin without it loads no particle code and refuses a `.star` file. The breaking-change note is written.
3. P2, the lazy `ParticlesFeature` with its stub's `owns`. Tests: a `.star` dropped on a Prism without the chunk loaded loads it, then the file. A session with particles does the same.
4. P1, node kinds, with tests.
5. P3 and P4: the inspector and the row menu.
6. P5: the unit cell row.
7. P6: targets, the API, then the inspector.
8. P7: loading choices.
9. P8: Simularium frames.
10. P9: the stories, `packages/prism/viewer/src/views/particles/_stories/particles.stories.tsx`:
   - each format;
   - the tree, with a list, its representations and its unit cell;
   - the inspector;
   - targets with 1crn;
   - a Simularium playback.
11. Write-up: the specs (`nodes.md`, `viewer-spec.md`, `parity.md`), the TODO, a `.v6/changelogs/` entry and `CHANGELOG.md`. Then this folder is deleted.

## Decisions (David, 2026-10-21)

- **Q1.** All of P0–P9 before clipping, with stories.
- **Q2.** Download, trim and commit the fixtures, with their licenses.
- **Q3.** Prism only. The classic `apps/viewer` drops particles (it does not compose `ParticlesRegistry`). It can opt in later if someone asks. This is recorded as a break.
- **Q4.** mmCIF particles are claimed by a content sniff in the particles feature's stub.
