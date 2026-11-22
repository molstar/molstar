# Plan: Web workers (the first in Mol\*)

## Status

- **Proposed**: no worker in the project yet. The first use: a dropped particle pack read off the main thread
  ([particle-tiles.md](particle-tiles.md) 10f2). What any worker needs (how it is built and found, its fallback, its
  messages, its tests) settled here once, before that use.
- **A minimal inline worker tried** (a dropped file's reads alone: its `Blob.slice` and `arrayBuffer` in the worker,
  the merging on the main thread, a fallback where none starts): the reads gone from the main thread (its idle time in
  the load's 4 s 270 → 900 ms) but the load's slow frames as before (p95 29-32 ms, long tasks about 200 ms): their
  cost the steps' (members filled for new targets, the merged object's passes after each, the uploads), not the
  reads. Not taken: a worker helps a phone's budget (the reads four times as long) more than a desktop's frames.

## Why

- **Measured**: the yeast's first tiles from a dropped pack (the tiles benchmark's still segment outside, at a
  display's pace): about 45 of 210 frames past a display's frame (about 29 ms), long tasks 140-280 ms in 4 s. Of
  that, the file's reading on the main thread: `Blob.slice` about 120 ms (a read's call), `Blob.arrayBuffer` about
  190 ms, their garbage about 100 ms; the steps' own stages about 8 ms of such a frame. Reads merged wider or kept
  as views do not change the frames (particle-tiles.md 10f2). A tile's decoding is views over its bytes: nothing.
- **Later uses** (not planned here): reads of a URL's ranges; parsing and other work that imports Mol\* code.

## Where a worker must work, or step aside

- **Builds**: the Viewer an IIFE (`molstar.js`, a classic script: no `import.meta`); Prism ES modules split into
  chunks (`scripts/esbuild/app.mjs`, `format: 'esm'`); the packages ES modules (`lib/`) bundled by their users' tools
  (Vite, webpack, esbuild, none).
- **Runtimes**: browsers; headless Node (`@molstar/plugin-headless`); tests (Vitest, Node, jsdom). No `Worker` in
  Node.
- **Pages**: an embedding page's content security policy may refuse a worker (`worker-src`, `blob:`).
- **Rule**: a worker is an optimization, never required: where it cannot start (none, refused, failed later) the
  same work on the main thread, the same results.

## How a worker is started: the options

1. **Inline**: the worker's code a function of its own (no imports, self-contained), its source a `Blob` URL
   (`new Worker(URL.createObjectURL(new Blob([...])))`). Works in every build and for every bundler (no entry, no
   asset); refused by a policy without `worker-src blob:`; its code cannot import Mol\* modules.
2. **Module worker by URL** (`new Worker(new URL('./x.worker.js', import.meta.url), { type: 'module' })`): its code
   imports Mol\* modules. Recognized by Vite and webpack 5, not by esbuild without a plugin; none in an IIFE build (no
   `import.meta`): the Viewer's would need its URL given.
3. **A URL given**: an app builds a worker entry and passes its URL (a plugin config); a library's user hosts it.

**Proposed**: 1 for work that needs no Mol\* code (the pack's reads: `Blob.slice`, `arrayBuffer`, the ranges' merging,
a few dozen lines), with the main thread as the fallback; 2 or 3 taken up when a worker needs Mol\* code, decided
then with the builds that need it.

## Design

- **The helper** (core, `util/worker.ts`): a worker from a self-contained function (`createInlineWorker`: undefined
  where `Worker` is none or its creation throws); requests and replies by id over a `MessagePort` (a worker's or a
  `MessageChannel`'s), their buffers transferred (no copy), an abort a message, an error a rejected reply; a worker
  failing (its `error` event: a policy refusing it after creation) rejects what is pending and tells its owner, who
  goes on without it. Ended (`terminate`) with its owner.
- **The pack's reads** (model `pack-stream.ts`): the ranges' merging as `batched` does it a pure function (shared by
  both: tested once); `ParticlePackSource.fromBlobInWorker(blob)`: the worker holds the `File` (passed by reference),
  merges the pending ranges, reads them on its thread, copies each range out and transfers its buffer; the main
  thread's source where it has no worker or loses it (pending reads asked again there). A source's `dispose` (new,
  optional) ends its worker, called as its tile set goes.
- **Where it is chosen** (plugin `particle-pack-source.ts`): a file's source the worker's, a URL's as now.

## Tests

- **Node**: the merging (ranges in, reads and parts out); the helper's requests, transfers, aborts and errors over a
  `MessageChannel` (Node has it, and structured cloning), no worker; the fallback where `Worker` is none.
- **A browser**: the tiles benchmark (the load's slow frames and long tasks, the reads in the worker); a Prism smoke
  run opening a pack by file.

## Steps

1. The helper and its tests (Node, `MessageChannel`).
2. The pack's reads in a worker, the fallback, the source's `dispose`; tests.
3. Measured by the tiles benchmark; the plans and specs (core's and model's `_spec` where they describe sources).
4. Later: a URL's reads (`fetch` in the worker); a module worker where Mol\* code is needed (option 2 or 3).

## Open questions (the maintainers')

- **The pattern**: inline workers for self-contained work, module workers when Mol\* code is needed: or one way for
  all from the start?
- **Policies**: do the embedders that matter (PDBe, RCSB, the Prism app's hosts) allow `worker-src blob:`? Without it
  inline workers never start there (the fallback: as today).
- **Node**: headless runs on the main thread, or `worker_threads` behind the same helper?
