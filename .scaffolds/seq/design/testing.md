# Testing strategy

## Intent

The author, 2026-11-17: "also what is our testing strategy?"

Tests must prove the promises the other designs make (the features, the edge cases, the performance targets, the
foundations' rules, the ports' contracts) and stay fast enough for the edit loop. The previous implementation's tests
were thorough and heavier than its code (about 119k lines against 103k, seven tiers); here the same promises are held
by fewer, more generic mechanisms.

## Decision

### What is tested, by what

| Promise | Held by |
|---|---|
| The language: validation, defaults, migrations, every issue code | per-level tests from the language description: each kind's examples valid, each documented misuse giving its code |
| Encodings and printing | properties: `decode(encode(x)) = x` for MessagePack and JSON, `json(msgpack(x)) = json(x)`, `parse(print(x)) = x`, over generated trees (typed columns of every kind, NaN, custom fields, extension kinds) |
| Patches | properties: `apply(a, diff(a, b)) = b`, `apply(apply(a, p), inverse(p)) = a`, over random edits; all-or-nothing on failure |
| Addresses | properties: `parse(print(a)) = a`; selection bounded; every numbering case (insertion codes, IMGT order, cyclic, reversed pieces); rects on every two-axis geometry |
| Maps | the same cases on every encoding of the same map (implicit, runs, expanded), equal answers; paths against a brute-force oracle on small graphs |
| Ports | one conformance suite per port, run on every adapter, built-in or an application's (§suites-next-to-ports), with fixtures, never the network |
| Layouts and the foundations' rules | the layout suite on every layout and preset: grid, one label column, aligned columns, density steps that don't move content, emphasis order, shapes |
| The interaction grammar | the grammar suite: each default binding and preset, on each layout and target kind, against the models in Node (no DOM) |
| Edge cases | each kept item a named test in its area's folder; the area's plan lists them, its record says each outcome ([edge-cases.md](edge-cases.md)) |
| Examples | every view: it validates, resolves and lays out with no unexpected issues and renders to SVG (quick); it opens in Storybook (medium); about six representative views keep a committed draw-list record and printed spec, compared as text |
| Rendering | the draw lists through the recording backend (quick); the Canvas2D and SVG backends' pixels compared (heavy) |
| The resolved view is general | the test-only second consumer draws its subset from the resolved view and the data store alone ([rendering.md](rendering.md)) |
| Performance | time bounds (`itWithin`), scaling ratios (cost follows pixels: ten times the data at most 1.2 times the frame; twice the pixels at most twice), allocation in the hot loops, heap growth on retire ([performance.md](performance.md)) |
| Size | a check on the built entries against the size targets |
| Parsers and readers | fuzzed (formats, addresses, CIGAR, the S-expression reader): never a throw, always issues |

### Fast by budget

The author: "we should keep fast test in mind too to implement this efficiently." The previous implementation's tests
slowed the edit loop; here speed is a requirement with numbers:

- **The quick lane for all of `packages/seq` runs in under 30 s** on the reference machine, each folder's tests in
  under 2 s, a test file in under 1 s (the review's budgets, accepted). A quick test over 100 ms is looked at. A test past its budget is fixed (waste) or moved with a quick representative
  (nature), in the change that made it slow.
- **No warm-ups, no loops of thousands, no wall clocks in the quick lane**; fixtures computed once per file, never per
  test; one viewer per file where a viewer is needed; property tests at a small, fixed number of seeded cases (more in
  the heavy lane).
- **Records stay few**: about six representative views keep a draw-list record, at one size, compared as text; the
  rest are checked by properties (no issues, bounded draw lists, hits where marks are), not by records.
- **The opt-in tiers stay opt-in**: allocation and heap guards run by hand at checkpoints, never in the edit loop or on
  every push.

### Where each runs: the repository's lanes

- **Quick** (every push, the edit loop; the whole of `packages/seq` in seconds): unit and property tests at small
  sizes, the conformance suites, the edge cases, every example's validate, resolve, layout and SVG render, the few records, the models
  and the grammar in Node, the layout suite on the default preset, short fuzz runs.
- **Medium** (every push, its own job): the layout suite on every preset, every example in both encodings and every
  preset, the Storybook smoke (every story opens without an error), the chrome's accessibility (axe), larger fixtures.
- **Heavy** (nightly, by hand, labelled PRs): time bounds and scaling ratios, full fuzz, pixel comparisons in Chrome,
  the render CLI in a process, the scale fixtures (a whole genome by region, 10M lazy rows).
- **Opt-in tiers, by hand at checkpoints** (the shared `test-harness`): allocation (`quickCase`) for the draw and hover
  loops only, heap guards on viewer retire. Not the previous implementation's seven tiers.

### How tests are written

- **Models in Node**: every model is tested without a DOM; views opt into jsdom per file; timers and frames through
  `DomHost`, patched in tests, so no test waits on a clock.
- **Records, not screenshots, in the quick lane**: draw-list records and printed specs are committed text, reviewed as
  diffs; a change to one is a change to look at, not to accept blindly.
- **Fixtures are small and committed**; large ones are generated in code (seeded), never stored; the network is never
  reached (the asset port's test adapter refuses it).
- **The previous implementation's answers are guidelines, not fixtures** (the author, 2026-11-17: "let's take them as
  guidelines. i think the testing in the old package got really slow and messed up the dev process a lot"): where a
  behavior is kept, its old answers (maps, projections, numbering, layout positions) are read when writing the test and
  the code, and the test states the expected answer itself, small and local. No recorded oracle runs.
- **Proportionate** (§proportionate-tests): a test costs about what its subject does; a slow test by nature moves lane
  and keeps a quick representative; a slow test by waste is fixed.
- **A story that composes more than one consumer has a headless twin** (§stories-are-tests, §whole-composition-tested):
  a test next to the story composes it as the story does (the viewer, its sources, a canvas and an export, two views),
  in Node on a fake clock. Why: what breaks between consumers (one releasing what the other reads, a signal one waits
  on that the other's change never fires) cannot be seen by each view's own tests; one such twin found three bugs.
  **The canvas against its export after a move**: the export lays out afresh from the store, the canvas from what it kept, so
  a pan and a zoom of every side, the canvas's ribbon count against the export's, finds a cache keyed too coarsely
  (`across-moves.test.ts`, FX28: a read kept its rows of the window it was first built for).
- **Every asynchronous source has a synthetic twin with an injectable clock**: a source runtime lands with a synthetic
  source of the same protocol whose delays run on an injected clock, used by its stories and its tests
  ([sources.md](sources.md)). Why: with it every asynchronous bug reproduces in milliseconds with no real waiting, and a
  story and its test see the same delays.

## Consequences

- A feature is `done` ([features.md](features.md)) when its view passes the three ways and its edge cases have
  outcomes; a performance target holds only once its test exists.
- Each port, layout and language ships its suite beside it, so an extension is tested by the same suite as the
  built-ins.

## In the specs

| Decision | Where |
|---|---|
| What is tested by what (language, encodings, patches, addresses, maps, ports, layouts, grammar, rendering, size, parsers) | design-only: the test plan; the tests are the record and each area's spec states its contract |
| Conformance suites per port and adapter, never the network | core › Protocol |
| Every example validates, resolves, lays out and renders; a few keep a record | examples › Tests |
| Performance targets: time bounds, scaling ratios, allocation, heap | design-only: the performance targets and the root AGENTS.md's tiers hold them |
| The quick lane's budgets (30 s the package, 2 s a folder, 1 s a file), no warm-ups or wall clocks | design-only: the package guide's Tests and the root AGENTS.md's lanes |
| Records few, at one size, compared as text, reviewed as diffs | examples › Tests |
| The opt-in tiers stay opt-in | design-only: the root AGENTS.md's "Not lanes" |
| Which lane each test goes to (quick, medium, heavy) | design-only: the root AGENTS.md's lanes |
| Models tested in Node, views opt into jsdom | viewer › Models |
| Fixtures small and committed, large ones generated from a seed | examples › Fixtures |
| The previous implementation's answers are guidelines, not fixtures | design-only: a working rule for writing tests |
| Proportionate tests: slow by nature moves lane, slow by waste is fixed | design-only: the root AGENTS.md's lanes |
| A story composing more than one consumer has a headless twin | examples › Tests |
| Every asynchronous source has a synthetic twin with an injectable clock | core › Asynchronous sources |
