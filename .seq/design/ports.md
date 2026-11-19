# Ports: functionality wired at composition

## Intent

The author, 2026-11-17: on archives, "i think we should support out of the box .. but via a "port" .. in general .. a
lot of the functionality should be wired up via "ports""; on assets, "reuse asset manger from
packages/core/src/util/assets.ts?"; after the review, its cuts accepted ("rest as recommended").

The core and the viewer reach nothing outside directly. What is one per viewer and touches the outside is a port the
application binds in the viewer's spec ([viewer-composition.md](viewer-composition.md)); what is a named provider (a
format, a map runtime, a backend) is a registry entry. Each has a default, so the viewer works out of the box.

## Decision

- **The rule**: one per viewer and outside-facing → a port (a factory with ownership in `spec.ports`); several, looked
  up by name → a registry entry ([extensibility.md](extensibility.md)). Nothing is both.
- **The ports of this version**:

  | Port | What | Default adapter |
  |---|---|---|
  | assets | read bytes by location | `@molstar/core/util/assets`'s `AssetManager`, wrapped by a policy (allowed hosts, size limits, no credentials by default, [security.md](security.md)); a host passes its own manager (a Mol\* plugin's), so downloads and files are shared, or a bare `{ reader }` (a test or a CLI has a reader, not a manager) |
  | text metrics | measure text for layout | `@molstar/core/util/text-measure` with a canvas in the browser, a bundled metrics table in Node ([rendering.md](rendering.md)) |
  | structure | bind structure systems to a structure viewer (hover and selection both ways, colors) | none; the Mol\* adapter in the showcases |
  | hover card | show the hovered place's card (`show(card, place)`, `hide()`): the card's content as data (a title, rows of a label and a value, the announcer's text), its place as boxes (the anchor to stay clear of, the pointer, the view, a suggested side from `placeCard`); a host draws it with its own components and styles (the author, 2026-11-19: "implement the hover card via a port too … so apps can easily provide their own styles") | `@molstar/ui`'s card (`hoverCardOverlay()`), passed by the default spec as the `hoverCard` behavior's `port` param |

- **A UI port a behavior owns lives on that behavior** (the author, 2026-11-19, HC): the hover card's port is the
  `hoverCard` behavior's `port` param, not `spec.ports`, so a viewer without the behavior carries no card; the default
  spec supplies the default adapter, and `disposeWithViewer: false` shares a host's adapter between viewers.
- **Registry entries, not ports**: format readers and writers (FASTA, alignments, tables, Newick, bigWig, SVG, PNG,
  FASTA export), map runtimes, backends (Canvas2D, SVG, recording), services (UniProt, PDBe, AlphaFold DB, in the
  showcases).
- **Not ports**:
  - **time, frames, timers**: `@molstar/ui/dom`'s `DomHost`, a module-level host the UI runtime reads directly; tests
    patch it as `@molstar/ui`'s do. A per-viewer host would diverge from the UI's own scheduling.
  - **archives**: supported out of the box as functions over the assets port (a spec with its data files and an
    index, read and written as zip with `@molstar/core/util/zip`, its entries resolved as assets), with limits on
    entries, sizes and paths (EC-91). The author asked for archives "via a port": the assets port is that port.
  - **clipboard, storage, links**: commands and plain functions over the platform, replaceable by binding other
    commands; they become ports only when a host needs a second implementation.
- **Later, when a profile or a host asks**: compute (a worker pool for summaries, alignment, parsing), with its
  boundary decided then (what crosses, columns transferred, the asset manager in a worker).
- **Suites only where a second adapter exists**: the assets port (the default and a test adapter that refuses the
  network), text metrics (canvas and table), and the format readers' shared contract.
- **The asset manager is reused, not rewritten**: it holds downloads and files ref-counted, shares one download across
  readers, keeps bytes once and inflates gzip once. What it lacks is range reads for indexed formats (bigWig): an
  extension of the manager (a page cache by asset and page, ref-counted, `Blob.slice` for files), a public API change
  proposed in [base-library.md](base-library.md).

## Consequences

- The core's dependencies stay `@molstar/core` plus its allowlist ([packages.md](packages.md)); an adapter that needs a
  platform or a library lives with the adapter.
- The Storybook's Extensibility page lists the ports and the registry kinds with their defaults.
