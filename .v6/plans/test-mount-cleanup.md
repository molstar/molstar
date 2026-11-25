# Mounted UI disposed between tests

Status: planned.

## Why

Test files share their worker (`isolate: false`, AGENTS.md, Tests), and jsdom's document with it. A view a test mounts
and never disposes keeps its scope, and with it its listeners on `document` and `window`: a later file's events reach it
and it wakes up (renders, attaches listeners, sets state). Seen twice, both as flakes in full runs only:

- `measurements-add.test.tsx` (a drag's lost releases) counted 44 listeners its drags did not add: shells left by
  `cards.test.tsx`, re-rendering on its pointer events (fixed in that file).
- The waste audit over budget: a real bug (a lazy tab's late load), but a leftover could cause the same.

Measured (the quick lane of `packages/prism`, `extensions`, `sequences`, `ai`, `apps`): 64 files leave 232 roots
mounted when they end (Prism UI 35 files, the viewer 21, the assistant 3, extensions 5; Seq\* none). Most mount small
views that never listen on the document, so they are harmless today; the next one that does is a flake in another file.

## Design: a registry of roots, disposed by the shared setup

- **Roots registered where the setup asks**: a root scope (no parent: `render` without one, `new Scope()`) registers
  itself in a set on `globalThis` under a symbol, when the set is there, and leaves it when disposed. The setup puts the
  set there; apps never do (a lookup per root, no import of the setup, nothing kept). In `Scope`'s constructor, not in
  `render`: tests make root scopes directly too (9 files: `new Scope()`).
- **Disposed by `vitest.setup.ts`**: the roots made during a test after it (`afterEach`), those made outside a test
  (`beforeAll`, the module) after the file (`afterAll`). The setup's `beforeEach` marks the running test, so a root made
  in a file's `beforeEach` is that test's. A test that disposes its own leaves nothing; disposing twice is not possible
  (the registry holds only live roots). A file that keeps a view across its tests mounts it in `beforeAll`, as the
  waste audit does, and keeps it to its end.
- **Cleanup, not failure**: the setup disposes what is left, as `@testing-library`'s auto-cleanup; it does not fail the
  file. Why: the leak is harmless to the file that leaves it and fatal only to a later one; disposing it removes the
  hazard wherever it is, and 64 files need no change. A report mode (an environment variable) lists the files that
  leave roots, for whoever wants them clean.
- **The DOM too**: the setup empties `document.body` after each file (jsdom only), as most files already do after each
  test: nodes left there are found by later files' queries.
- **Both runtimes**: Prism's (`packages/prism/ui/src/dom`) and `@molstar/ui/dom` (`packages/ui`), the same registry
  symbol, so the setup knows neither package.

Rejected: a lint check that a test's `render(…)` result is used. Tests mount through helpers (`mount`, `headlessViewer`),
keep views on purpose across a file, and make scopes without `render`; a runtime registry sees what is actually left.

## Steps

1. **The registry** in both runtimes' root scopes, and the setup's hooks (`beforeEach` marks the test; `afterEach` and
   `afterAll` dispose). Tested: a root made in a test disposed after it, one made in `beforeAll` kept to the file's end,
   one the test disposed not disposed again; apps without the set unchanged (`scope.ts`'s tests).
2. **A file of its own as the guard's test**: mounts a view listening on the document and leaves it; the next test sees
   its listener gone (the regression `cards.test.tsx` caused).
3. **The report mode** (`MOLSTAR_REPORT_ROOTS=1`): after each file, the files and mount sites of the roots the setup
   disposed; measured once against the 64 above.
4. **Docs**: AGENTS.md (Tests: what a file restores, and what the setup disposes for it), Prism's AGENTS.md (the
   timing-sensitive checks: a leftover is disposed, a failure under load is the file's own).
