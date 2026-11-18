# .seq: the development scaffold of Seq\* (packages/seq)

`packages/seq` is a new implementation of Seq\*: SeqViewSpec (a JSON language for sequence views over one coordinate
model), its model, its query language and its viewer, built on `@molstar/core/reactive` and `@molstar/ui`. This folder
is its agent context while it is being built: the principles it follows, the sources it is drawn from, the agreed
designs, the work folders and the records.

**It is temporary.** When the package stands on its own, what must outlive this folder (the rules the code follows,
the decisions, what the format means) is moved into the packages' `_spec/` and `.seq/` is deleted. So:

- nothing outside `.seq/` links into it except the one pointer in `packages/seq/AGENTS.md`;
- code, comments, tests and specs never cite it (no `.seq/` path, no `§` id of the principles in code);
- `packages/seq` is written as if `packages/sequences` (the previous implementation) did not exist: it may be read as
  a source here, never named in the package's code or specs, since it will be deleted.

The rest of the repository follows the root `AGENTS.md` and `.v6/`; both apply here too, and where they speak about the
same thing (package rules, tests, headers, changelog) the root wins.

## Where things are

| Path | What |
|---|---|
| [status.md](status.md) | **start here**: what is built, what is open, what is next, which design to read |
| [principles/README.md](principles/README.md) | the manifest: which modules apply, the folder map, the adaptations |
| `oz/` | the principles package, copied unedited ([its README](oz/README.md)) |
| `principles/{queries,agent-interface}/` | the project's own modules, from meso\* |
| [sources/](sources/README.md) | what Seq\* is drawn from, each reconciled against the principles |
| [design/](design/README.md) | agreed designs, one file per topic; if it is here it is agreed |
| [records/](records/README.md) | one permanent record per execution session |
| `execution/<work-id>/` | plan, loop, state and drafts of work in flight; committed and kept, never cited from outside |
| `_skills/seq/` | the instance skill |
| `_checks/` | this scaffold's own checks, beside the package's (Checks, below) |

## How work happens here

0. Read [status.md](status.md), then the [instance skill](_skills/seq/SKILL.md).
1. Read the `principles` skill ([oz/_skills/principles/SKILL.md](oz/_skills/principles/SKILL.md)), then the router
   ([oz/_skills/principles-router/SKILL.md](oz/_skills/principles-router/SKILL.md)) with the manifest for the modules
   of the area you touch. The skills are not linked into the harness (so deleting `.seq/` leaves nothing); read them as
   files.
2. Design first: a draft in `execution/<work-id>/`, agreed into `design/<topic>.md`; then plan, loop and `state.md` in
   the execution folder (read `state.md` first when resuming).
3. Code, spec and design disagree: write a drift note, don't bend either side. Changes to `principles/`, `sources/` or
   an agreed design are the author's call (David): bring the evidence and a recommendation.
4. Before calling work done, run the checks yourself, the repository's and these; a report, yours or a subagent's, is
   not proof.

Everything in `.seq/` is agent context: never auto-formatted (`.prettierignore`), and machine-invariant (no local
paths in committed files).

## Checks

They read the git index: stage first (`git add .seq packages/seq`).

```sh
node .seq/oz/_checks/principles-ids.ts --package .seq/oz --instance .seq --allow '^(?!\.seq/)'
node .seq/oz/_checks/paths-check.ts --allow '^(?!\.seq/|packages/seq/)'
node .seq/oz/_checks/merge-gate.ts --package .seq/oz --instance .seq
node .seq/_checks/ec-ids.ts
node .seq/_checks/closed-work.ts
```

The last two are this scaffold's own: every `EC-<n>` cited in `packages/seq` names an item of the edge-case catalog
([sources/seqstar/edge-cases.md](sources/seqstar/edge-cases.md); [design/edge-cases.md](design/edge-cases.md)); and a
closed work folder has a `foldback.md` the author approved ([principles](principles/README.md), "The fold-back is a step
of every close").

The `--allow` patterns limit the id and path checks to what this scaffold governs: the rest of the repository uses
the section sign for its own spec sections (`harness.md`, section 6).
