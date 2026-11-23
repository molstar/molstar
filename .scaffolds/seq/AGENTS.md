# .scaffolds/seq: the development scaffold of Seq\* (packages/seq)

`packages/seq` is a new implementation of Seq\*: SeqViewSpec (a JSON language for sequence views over one coordinate
model), its model, its query language and its viewer, built on `@molstar/core/reactive` and `@molstar/ui`. This folder
is its agent context while it is being built: the principles it follows, the sources it is drawn from, the agreed
designs, the work folders and the records.

**It is kept, and never cited by code.** What the code follows (the rules, the decisions, what the format means) is in
the packages' `_spec/` (each design's `## In the specs` table says where). This folder keeps the reasons, the history
and the work, and may still be edited, the designs especially (the author, 2026-11-23: "i think we can still edit if
needed"). So:

- nothing outside `.scaffolds/` links into it except the pointer in `packages/seq/AGENTS.md` and the root `AGENTS.md`'s
  table;
- code, comments, tests and specs never cite it (no `.scaffolds/seq/` path, no `§` id of the principles in code);
- `packages/seq` is written as if `packages/sequences` (the previous implementation) did not exist: it may be read as
  a source here, never named in the package's code or specs.

The rest of the repository follows the root `AGENTS.md` and `.v6/`; both apply here too, and where they speak about the
same thing (package rules, tests, headers, changelog) the root wins.

## Where things are

| Path | What |
|---|---|
| [status.md](status.md) | **start here**: what is built, what is open, what is next, which design to read |
| [principles/README.md](principles/README.md) | the manifest: which modules apply, the folder map, the adaptations |
| `../principles/` | the principles package, copied unedited ([its README](../principles/README.md)) |
| `principles/{queries,agent-interface}/` | the project's own modules, from meso\* |
| [sources/](sources/README.md) | what Seq\* is drawn from, each reconciled against the principles |
| [design/](design/README.md) | agreed designs, one file per topic; if it is here it is agreed |
| [records/](records/README.md) | one permanent record per execution session |
| `execution/<work-id>/` | plan, loop, state and drafts of work in flight; committed and kept, never cited from outside; stays after the close, never moved or deleted unless the author says so (the record alone goes to `records/`) |
| `_skills/seq/` | the instance skill |
| `_checks/` | this scaffold's own check (`ec-ids.ts`) and its parts of the method's checks: the spec map, the scopes, the baselines (Checks, below) |

## How work happens here

0. Read [status.md](status.md), then the [instance skill](_skills/seq/SKILL.md).
1. Read the `principles` skill ([../principles/_skills/principles/SKILL.md](../principles/_skills/principles/SKILL.md)), then the router
   ([../principles/_skills/principles-router/SKILL.md](../principles/_skills/principles-router/SKILL.md)) with the manifest for the modules
   of the area you touch. The skills are not linked into the harness (so deleting `.scaffolds/seq/` leaves nothing); read them as
   files.
2. Design first: a draft in `execution/<work-id>/`, agreed into `design/<topic>.md`; then plan, loop and `state.md` in
   the execution folder (read `state.md` first when resuming).
3. Code, spec and design disagree: write a drift note, don't bend either side. Changes to `principles/`, `sources/` or
   an agreed design are the author's call (David): bring the evidence and a recommendation.
4. Before calling work done, run the checks yourself, the repository's and these; a report, yours or a subagent's, is
   not proof.

Everything in `.scaffolds/seq/` is agent context: never auto-formatted (`.prettierignore`), and machine-invariant (no local
paths in committed files).

## Checks

They read the git index: stage first (`git add .scaffolds packages/seq`).

```sh
node .scaffolds/principles/_checks/principles-ids.ts --package .scaffolds/principles --instance .scaffolds/seq --allow '^(?!\.scaffolds/)'
node .scaffolds/principles/_checks/paths-check.ts --allow '^(?!\.scaffolds/|packages/seq/)'
node .scaffolds/principles/_checks/merge-gate.ts --package .scaffolds/principles --instance .scaffolds/seq
node .scaffolds/_method/_checks/closed-work.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/defaults-designed.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/designs-in-specs.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/no-design-citations.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/file-size.ts --instance .scaffolds/seq
node .scaffolds/seq/_checks/ec-ids.ts
```

The first three are the principles package's, the next five the method's ([../_method/README.md](../_method/README.md),
Checks), each reading its Seq\* part from a file here in `_checks/`: `closed-work.json` (the folders closed before the
fold-back rule), `defaults-designed.json` (the designs exempt, only shrinks), `designs-in-specs.json` (the packages'
specs by key), `no-design-citations.json` (Seq\*'s packages, the one pointer), `file-size.json` (`packages/seq`) with
`file-size-baseline.json` (only shrinks, lowered by `--update`). A loop also runs the method's `features-planned.ts
--instance .scaffolds/seq --plan <the slice's plan.md>`.

`ec-ids.ts` is this scaffold's own: every `EC-<n>` cited in `packages/seq` names an item of the edge-case catalog
([sources/seqstar/edge-cases.md](sources/seqstar/edge-cases.md); [design/edge-cases.md](design/edge-cases.md)), and
a feature marked `done` in [design/features.md](design/features.md) has each item it covers named by a test or given an
outcome (rows done before the rule in `_checks/ec-coverage-baseline.json`, which only shrinks; `--worklist` prints the
built rows' items no test names, the close's list).

The `--allow` patterns limit the id and path checks to what this scaffold governs: the rest of the repository uses
the section sign for its own spec sections (`harness.md`, section 6).
