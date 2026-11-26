# Mounted UI disposed between tests

Status: done.

## Why

Test files share their worker (`isolate: false`, AGENTS.md, Tests), and jsdom's document with it. A view a test mounts
and never disposes keeps its scope, and with it its listeners on `document` and `window`: a later file's events reach it
and it wakes up (renders, attaches listeners, sets state). Seen as a flake in full runs only:
`measurements-add.test.tsx` (a drag's lost releases) counted 44 listeners its drags did not add, from shells
`cards.test.tsx` left re-rendering on its pointer events.

Measured (the quick lane of `packages/prism`, `extensions`, `sequences`, `ai`, `apps`): 64 files leave 232 roots
mounted when they end (Prism UI 35 files, the viewer 21, the assistant 3, extensions 5; Seq\* none). Most mount small
views that never listen on the document, harmless today; the next one that does is a flake in another file.

## Design: a registry of roots, disposed by the shared setup

- **Roots registered where the setup asks**: a root scope (no parent: `render` without one, `new Scope()`) of either
  DOM runtime (Prism's `dom/scope.ts`, `@molstar/ui/dom`'s `ViewScope`) adds itself to a registry on `globalThis` under
  `Symbol.for('molstar.test.rootScopes')`, when there is one, and leaves it when disposed. The setup puts it there; apps
  never do (one lookup per root scope, nothing kept; Prism's allocation guards unchanged). In the scopes' constructors,
  not in `render`: tests make root scopes directly too. The setup imports neither runtime.
- **Disposed after each file** (`vitest.setup.ts`, `afterAll`, after the file's own hooks): the roots left, newest
  first, each on its own; then `document.body` emptied (jsdom). After the file, not each test: files keep views across
  their tests (a viewer made in `beforeAll`), and views make root scopes while a test uses them; disposing those after
  the test would break the next. Within a file, its tests always run together and in order.
- **Cleanup, not failure**, as `@testing-library`'s auto-cleanup: the leak is harmless to the file that leaves it and
  fatal only to a later one; disposing it removes the hazard wherever it is, and the 64 files need no change.
  `MOLSTAR_REPORT_ROOTS=1` writes where the roots left were made, to stderr (the console of a passing file is not
  shown).

Rejected: a lint check that a test's `render(…)` result is used. Tests mount through helpers (`mount`,
`headlessViewer`), keep views on purpose across a file, and make scopes without `render`; a runtime registry sees what
is actually left.

## Steps

1. **The registry and the setup** (done). Tested (`root-scopes.test.ts` in both runtimes): a root registered while
   live, its children not, gone when disposed, disposing again changing nothing; no registry (apps) keeping nothing;
   both fail without the registration. All quick-lane files of the packages above pass.
2. **Across files** (done, by hand): `cards.test.tsx` as before its fix (its shells left mounted) then
   `measurements-add.test.tsx` in one worker, the two orders that failed: both pass.
3. **The report mode** (done): `render.test.tsx` and `control-flow.test.tsx` report the 16 and 12 roots measured above.
4. **Docs** (done): AGENTS.md (Tests: what the setup disposes for a file), Prism's AGENTS.md.
