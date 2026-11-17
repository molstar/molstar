# seq: the principles here

`packages/seq` is a new implementation of Seq\*: SeqViewSpec, the language for sequence views, with its model and its
viewer, in the Mol\* v6 repository. `.seq/` is its instance of the oz principles package, and this README is the
manifest the `principles-router` skill reads. It governs `packages/seq/` and `.seq/` only; the rest of the repository
follows its root `AGENTS.md` and `.v6/`, which win where both speak about the same thing.

**`.seq/` is a development scaffold**, deleted once the package stands on its own. Nothing outside `.seq/` links into
it except the one pointer in `packages/seq/AGENTS.md`, and code never cites it. What must outlive it (the rules the
code follows, the decisions, what the format means) moves into the packages' `_spec/` folders before then.

## Provenance

- The package is `.seq/oz/`: the oz principles (`.principles/` in the oz exploration repository), copied at commit
  `7d03f60` (2026-11-17), unedited. An update from it is a merge from that base.
- `queries/` and `agent-interface/` here: from the meso\* instance of the same package (`.meso/principles/` in the
  molstar-ludo repository), copied at commit `84c41cd84` (2026-11-17). Their links were repointed at `.seq/oz/`, and
  queries' rule 14 (§bounded-work) became a pointer, the package's representations module now holding it.

## Folder map

| The package says | Here |
|---|---|
| `.principles/` | `.seq/oz/` |
| `.design/`, `.sources/`, `.records/` | `.seq/design/`, `.seq/sources/`, `.seq/records/` |
| `.execution/`, `.tmp/` | `.seq/execution/` (committed and kept, never cited from outside), `.seq/tmp/` (gitignored) |
| `_spec/` | each package's own `_spec/`, as everywhere in Mol\* |
| skills linked from `.agents/skills/` | not linked: the skill files are read directly ([.seq/AGENTS.md](../AGENTS.md)), so removing `.seq/` leaves nothing behind |
| the checks: id, path, merge | `node .seq/oz/_checks/<check>.ts --package .seq/oz --instance .seq` ([.seq/AGENTS.md](../AGENTS.md), Checks) |

## Modules

The router picks from the "Read before" column. A module of this project's own named like a package module extends
it: read the package's first.

### From the package

| Module | Kind | Read before | Covers |
|---|---|---|---|
| [_core](../oz/_core/README.md) | rule set | any non-trivial work (the `principles` skill loads it) | how the principles work and how work is run |
| [architecture](../oz/architecture/README.md) | rule set | packages, ports, wiring, errors, commands, performance | ports, API surfaces (sync, async, streaming), composition root, error codes, commands, performance |
| [representations](../oz/representations/README.md) | rule set | the language, documents, ids, patches, derived forms | declarative state as recursive, streamable trees; committed vs transient state; addresses, hashing; derived forms; encodings |
| [ui](../oz/ui/README.md) | rule set | the viewer's models and views, Storybook | frontend state and logic, design system and Storybook, views as representations |
| [testing](../oz/testing/README.md) | rule set | tests, conformance suites | test tiers (here the repository's lanes), conformance library |

### The project's own

| Module | Kind | Read before | Covers |
|---|---|---|---|
| [queries](queries/README.md) | rule set | the query language, locators, coordinate spaces, mappings | query DSLs, the semantic layer (concepts, spaces, mappings), embedded languages (from meso\*) |
| [agent-interface](agent-interface/README.md) | rule set | exposing the language to agents | skills, the `<domain>_read`/`_write` tool surface taking expressions, outline-first reads, affordances in results (from meso\*) |

## Adaptations

- **workspace**, **layout**: not adopted. The repository's toolchain, tree and package rules (root `AGENTS.md`,
  `check:workspace`) win.
- **architecture**: its layout (`lib/`, `packages/<product>/*`, `apps/`) maps onto Mol\*'s. The libraries
  `packages/{core,io,model,graphics}` and `packages/ui` play `lib/`, and Seq\* is the product (`packages/seq/*`). The
  root package rules win where they differ. Kept for ports, commands, errors and performance.
- **testing**: its tiers are the repository's three lanes (quick, medium, heavy; root `AGENTS.md`). Allocation and
  heap guards only for the draw and hover hot loops, once they exist. Conformance suites (§conformance-suites) are
  kept.
- **ui**: the runtime and design system are `@molstar/core/reactive` and `@molstar/ui` (no rxjs; React only as an
  optional adapter); §frontend-state is taken as written.
- **representations**: one tree, encoded as MessagePack (the default) or JSON, with an S-expression printer (and
  reader) for inspection, never a file type of its own (the author, 2026-11-17). This is §one-recursive-form's "the spelling is the
  project's choice", made here.
- **_core, `.execution/`**: committed and kept after the close, as in meso\*; §execution-unreferenced stands.
- **_core, followups**: open items go to `packages/seq/TODO.md` (the root `AGENTS.md` sends them to a package's
  TODO); an entry says its context, why it was deferred and what done looks like.
- **Headers**: every source file carries the Mol\* header with the year and `@author` lines only (the author,
  2026-11-17); otherwise the root `AGENTS.md` header rule.
- **Bootstrap patterns**: none taken in yet. The reactive-models pattern is replaced by `@molstar/core/reactive`
  (`ReactiveModel`, `State`, `derive`), which the host already provides.

## Open items

`packages/seq/TODO.md` once the package exists; until then [.seq/status.md](../status.md).
