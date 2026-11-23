# seq: the principles here

`packages/seq` is a new implementation of Seq\*: SeqViewSpec, the language for sequence views, with its model and its
viewer, in the Mol\* v6 repository. `.scaffolds/seq/` is its instance of the principles package, and this README is the
manifest the `principles-router` skill reads. It governs `packages/seq/` and `.scaffolds/seq/` only; the rest of the repository
follows its root `AGENTS.md` and `.v6/`, which win where both speak about the same thing.

**`.scaffolds/seq/` is a development scaffold**, kept, and never cited by code ([../AGENTS.md](../AGENTS.md)); what
the code follows is in the packages' `_spec/`.

## Provenance

- The package is `.scaffolds/principles/`: the principles (`.principles/` in the upstream exploration repository), copied at commit
  `7d03f60` (2026-11-17), re-synced at `c1ccd14` (2026-11-23: its rule budget, a lighter `_core`, Seq\*'s slices 3 to 5
  folded back). An update from it is a merge from that base. Edited on each sync: its `workspace` and `layout`
  modules removed (not adopted here) and their citations reworded, at the review's recommendation.
- `queries/` and `agent-interface/` here: from the meso\* instance of the same package (`.meso/principles/` in the
  molstar-ludo repository), copied at commit `84c41cd84` (2026-11-17). Their links were repointed at `.scaffolds/principles/`, and
  queries' rule 14 (§bounded-work) became a pointer, the package's representations module now holding it.

## Folder map

| The package says | Here |
|---|---|
| `.principles/` | `.scaffolds/principles/` |
| `.design/`, `.sources/`, `.records/` | `.scaffolds/seq/design/`, `.scaffolds/seq/sources/`, `.scaffolds/seq/records/` |
| `.execution/`, `.tmp/` | `.scaffolds/seq/execution/` (committed and kept, never cited from outside), `.scaffolds/seq/tmp/` (gitignored) |
| `_spec/` | each package's own `_spec/`, as everywhere in Mol\* |
| skills linked from `.agents/skills/` | not linked: the skill files are read directly ([.scaffolds/seq/AGENTS.md](../AGENTS.md)), so removing `.scaffolds/seq/` leaves nothing behind |
| the checks: id, path, merge | `node .scaffolds/principles/_checks/<check>.ts --package .scaffolds/principles --instance .scaffolds/seq` ([.scaffolds/seq/AGENTS.md](../AGENTS.md), Checks) |

## Modules

The router picks from the "Read before" column. A module of this project's own named like a package module extends
it: read the package's first.

### From the package

| Module | Kind | Read before | Covers |
|---|---|---|---|
| [_core](../../principles/_core/README.md) | rule set | any non-trivial work (the `principles` skill loads it) | how the principles work and how work is run |
| [architecture](../../principles/architecture/README.md) | rule set | packages, ports, wiring, errors, commands, performance | ports, API surfaces (sync, async, streaming), composition root, error codes, commands, performance |
| [representations](../../principles/representations/README.md) | rule set | the language, documents, ids, patches, derived forms | declarative state as recursive, streamable trees; committed vs transient state; addresses, hashing; derived forms; encodings |
| [ui](../../principles/ui/README.md) | rule set | the viewer's models and views, Storybook | frontend state and logic, design system and Storybook, views as representations |
| [testing](../../principles/testing/README.md) | rule set | tests, conformance suites | test tiers (here the repository's lanes), conformance library |

### The project's own

| Module | Kind | Read before | Covers |
|---|---|---|---|
| [queries](queries/README.md) | rule set | the query language, locators, coordinate spaces, mappings | query DSLs, the semantic layer (concepts, spaces, mappings), embedded languages (from meso\*) |
| [agent-interface](agent-interface/README.md) | rule set | exposing the language to agents | skills, the `<domain>_read`/`_write` tool surface taking expressions, outline-first reads, affordances in results (from meso\*) |

## Adaptations

What is Seq\*'s and in force now. The package's rules hold as written except where a line here says how Seq\*
applies them; the evidence behind each line is in the [records](../records/README.md).

- **workspace**, **layout**: not adopted, and removed from `.scaffolds/principles/`; citations of their ids read "not
  adopted here". The root `AGENTS.md` and `check:workspace` govern the toolchain, tree and package rules.
- **architecture**: `packages/{core,io,model,graphics}` and `packages/ui` play `lib/`; Seq\* is the product
  (`packages/seq/*`). The root package rules win where they differ.
- **testing**: the tiers are the repository's three lanes (quick, medium, heavy; root `AGENTS.md`); allocation and heap
  guards only on the draw and hover hot loops.
- **The heavy lane runs alone, at two workers** {#heavy-lane-alone}: no other test run beside it, `maxWorkers: 2` in
  heavy mode (`vitest.config.ts`), and a bound under 50 ms runs alone in its worker and is never raised to pass under
  contention (§test-tiers). Why: bounds failed in the full heavy lane under load and passed alone.
- **ui**: the runtime and design system are `@molstar/core/reactive` and `@molstar/ui` (no rxjs; React only as an
  optional adapter). No bootstrap pattern is taken in: `@molstar/core/reactive` (`ReactiveModel`, `State`, `derive`)
  replaces the reactive-models pattern.
- **representations**: one tree, JSON by default or MessagePack, with an S-expression printer and reader for
  inspection, never a file type of its own (§one-recursive-form's spelling, chosen here).
- **Followups** (§followups) go to `packages/seq/TODO.md`: context, why deferred, what done looks like.
- **Headers**: the Mol\* header with the year and `@author` lines only; otherwise the root `AGENTS.md` rule.
- **Lane models** {#model-by-scope}: §model-per-task, with Sonnet as the fast model.
- **Configuration tables** {#config-tables-checked}: each package's `_spec` keeps one Configuration table (item,
  default, why, alternatives), tested both ways against the default configuration (`src/_test/config-table.test.ts`
  in `seqx-core` and `seqx-viewer`); a design's Configuration row names an item the code has, or says "not built"
  (§mechanical-drift-checks). Why: copied by hand, rows went missing and designs named items never built.
- **Edge cases by name** {#ec-tested-by-name}: a test of a catalog item names its `EC-<n>` (§honest-tests);
  `.scaffolds/seq/_checks/ec-ids.ts` resolves every id cited in `packages/seq` against
  [sources/seqstar/edge-cases.md](../sources/seqstar/edge-cases.md), fails a `done` feature whose items have neither a
  test nor an outcome, and `--worklist` lists a built row's untested items (§mechanical-drift-checks). Why: the ids are
  bare numbers, so code cites them without a scaffold path.
- **The Features list moves with the merges** {#features-move-with-merges}: the lane that builds a feature moves its
  row in [design/features.md](../design/features.md) and its line on the Storybook's Design page "Features" (with its
  story's link) in its own change (§mechanical-drift-checks).
- **The Storybook's Design pages follow the designs**: a design change a reader of the package needs (a new design, a
  rule, a decision, a table) is folded into the matching page under `Seq* (next)/Design` in the same change, or a page
  is added; the pages never cite `.scaffolds/seq/`. Why: they are what readers see of the designs (the author,
  2026-11-18: "we should be folding the things we are adding to the designs section in the storybook too").
- **Interaction matrix** (§interaction-matrix): kept in [design/ux-foundations.md](../design/ux-foundations.md).
- **Defaults designed** (§defaults-designed): checked by `.scaffolds/_method/_checks/defaults-designed.ts`, its exempt
  list in `_checks/defaults-designed.json`.
- **Mocks** {#mocks-at-design}: §interactive-mocks, with the mocks in the Storybook under `Seq* (next)/Design/Mocks/`
  (not a `prototypes/` package), built by a fast lane beside the draft; no package imports a mock, and the merge that
  lands the built view deletes it.
- **A fold-back at every close** (tightens §reconcile-before-merge, which allows "Fold-back: none"): every close writes
  `foldback.md`, the author's decision on its own line (the method's template; `closed-work.ts` checks it).
- **The novice pass is standing** (tightens §used-before-done's "may"): it runs in every slice that draws, after the
  last drawing lane, and the orchestrator confirms its bugs before its findings are placed.

Tombstones (closed work folders cite them). Moved within the instance:

- {#heavy-bounds-isolated} *moved to §heavy-lane-alone*

Moved to the package:

- {#drift-notes-are-gaps} *moved to §gaps-queue*
- {#seeded-generators-checked} *moved to §suites-next-to-ports*
- {#test-expectation-explained} *moved to §honest-tests*
- {#matrix-cell-by-test} *moved to §interaction-matrix*
- {#use-continuous} *moved to §used-before-done*
- {#hot-loop-representatives} *moved to §performance-tested*
- {#findings-are-classes} *moved to §drift-notes*
- {#engine-not-game} *moved to §application-vs-base*

## For the package's home

- §honest-tests, optional at gates: a script that reverts a fix's non-test files on a scratch worktree and expects
  its new test to fail.

## Open items

`packages/seq/TODO.md` (§followups).
