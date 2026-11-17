# Ports: functionality wired at composition

## Intent

The author, 2026-11-17: on archives, "i think we should support out of the box .. but via a "port" .. in general .. a
lot of the functionality should be wired up via "ports"".

The core and the viewer do no I/O and reach no platform directly. Everything that touches the outside (files,
network, archives, fonts, time, storage, the clipboard, a structure viewer, threads) goes through a port: a typed
interface the application binds when it composes the viewer (§compose-through-ports, §composition-root). Most ports
ship a default adapter, so the viewer works out of the box; any of them can be replaced.

## Decision

- **A port is an interface plus a conformance suite** (§suites-next-to-ports): every adapter, built-in or an
  application's, passes the same suite, run against fixtures without the network.
- **Defaults out of the box**: `createViewer()` with no arguments binds the default adapters; an application passes
  its own to replace any of them. A port with no adapter makes what needs it unavailable, with an issue, never a crash.
- **Ports are bound in the viewer's spec** (`spec.ports`, factories with ownership, [viewer-composition.md](viewer-composition.md)),
  so what a viewer can reach is visible in one place; runtimes that are providers by name (formats, maps, services) are
  registry entries instead.
- **The ports, first list**:

  | Port | What it does | Default adapter |
  |---|---|---|
  | assets (fetch) | read bytes by location and range, with a policy (allowed hosts, size limits, credentials) | `@molstar/core/util/assets`'s `AssetManager`, extended with range reads; a host passes its own (a Mol\* plugin's), so downloads and files are shared |
  | archive | read and write a bundle: a spec with its data files, resolved by relative location | a zip archive (read and write), a directory, a map in memory |
  | format runtimes | a data format onto the one source protocol ([sources.md](sources.md)) | FASTA, alignments, tables, Newick, MessagePack and JSON; bigWig as the proof of a genome format |
  | map runtimes | compute a map (align, a SIFTS lookup, a gene model to codons) ([maps.md](maps.md)) | alignment of two sequences; others an application adds |
  | structure runtime | bind structure systems to a structure viewer, link hover and selection, color by tracks | none in the core; the Mol\* adapter in the showcases |
  | text metrics | measure text for layout | a canvas in the browser, a bundled metrics table in Node ([rendering.md](rendering.md)) |
  | backends | draw draw lists | Canvas2D, SVG, recording |
  | writers | export (a selection as FASTA, a view as SVG or PNG, a spec as an archive) | SVG, PNG, FASTA, archive |
  | host | time, frames, timers, idle callbacks (so models run in Node and tests control time) | `@molstar/ui/dom`'s host |
  | compute | run heavy work off the main thread (summaries, alignment, parsing) with cancellation | a worker pool in the browser, inline in Node |
  | storage | keep user settings and drafts | in memory; local storage in the browser |
  | clipboard | copy links, addresses, sequences | the browser's clipboard; in memory in Node |
  | links | open an address from a link, hand an address to the host (a URL, a markdown renderer) | Mol\*'s link grammar ([addressing.md](addressing.md)) |
  | services | named data services (UniProt, PDBe, AlphaFold DB) as source runtimes | none in the core; the showcases compose them |

- **Archives are supported out of the box**, through the archive port: a spec (MessagePack, or JSON) with its data
  files and a small index, sources inside it named by relative location; read and written by the default zip adapter;
  limits on entries, sizes and paths enforced by the port (EC-91), never by the caller; an archive's entries resolve as assets of the same manager.
- **The asset manager is reused, not rewritten** (the author, 2026-11-17: "reuse asset manger from
  packages/core/src/util/assets.ts?"): `AssetManager` already holds downloads and files ref-counted, shares one download
  across readers and response types, keeps bytes once, inflates gzip once, and is what a Mol\* plugin uses. Seq\*'s
  sources resolve through it, so a structure file the Mol\* showcase opens is downloaded once for both. What it lacks
  is range reads for indexed formats: an extension of the manager (a page cache by asset and page, ref-counted like its
  entries, `Blob.slice` for files), proposed in [base-library.md](base-library.md).
- **Reaching the outside world is a policy**: the asset port's adapter decides which hosts and sizes are allowed;
  tests bind an adapter that refuses the network.

## Consequences

- The core's dependencies stay `@molstar/core` only; every adapter that needs a platform or a library (a zip
  implementation, a canvas, Mol\*) lives with the adapter, loaded only when bound.
- The Storybook's Extensibility page lists the ports with their defaults; each port's suite runs in the quick lane.
- The size targets ([performance.md](performance.md)) count the core without adapters; each adapter has its own.
