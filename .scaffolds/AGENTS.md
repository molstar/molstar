# .scaffolds: development scaffolds

A **development scaffold** is the agent context a piece of the repository is built in: the principles it follows,
the sources it is drawn from, its agreed designs, its work folders and its records. It stands beside the code while
the code is built and is kept after, as the reasons and the history behind it. What the code follows (its rules, its
decisions, what its formats mean) moves into the packages' `_spec/`; the scaffold keeps the why and the how.

The rest of the repository follows the root `AGENTS.md` and `.v6/`; they apply inside a scaffold too, and where both
speak about the same thing (package rules, tests, headers, the changelog) the root wins.

## Layout

| Path | What |
|---|---|
| [principles/](principles/README.md) | the principles package, shared by every scaffold: rule modules, patterns, skills, the three package checks; copied from its upstream, edited only as its manifest's provenance says |
| [_method/](_method/README.md) | the method: templates for a work folder's files and a design, and the checks every scaffold runs on itself (`--instance .scaffolds/<name>`) |
| [seq/](seq/AGENTS.md) | Seq\*'s scaffold (`packages/seq`): its manifest and own modules, sources, designs, records, work folders and its own check |

One folder per scaffold, `.scaffolds/<name>/`, laid out as the method's README says: `AGENTS.md`, `status.md`,
`principles/README.md` (the manifest), `_skills/<name>/`, `sources/`, `design/`, `records/`, `execution/`, `_checks/`.

## A scaffold's life

- **Started** from `_method/` ([its README](_method/README.md), "A scaffold's life"): the instance folder, its manifest
  from the package's template, its checks' instance files, a pointer in the root `AGENTS.md`'s table and in the
  governed package's `AGENTS.md`.
- **Run** as slices: designs agreed first, then a work folder per slice from the method's plan, loop and state, one
  brief per lane.
- **Closed** by the loop's Finalize: reconcile, the fold-back with the author's decision, the record, the next plan.
  A work folder stays where it is after its close; only the record goes to `records/`.
- **Kept** when the governed packages stand on their own (the author, 2026-11-23: "originally i planned to remove it,
  but i think we should reconcile it instead and put it in some common place, where next similar "scaffolds" will
  live").

## Rules

- **Kept and editable.** A closed scaffold is not deleted and not frozen: its files may still be edited, the designs
  especially (the author, 2026-11-23: "i think we can still edit if needed (and prolly will, especialyl some of the
  designs)"). A closed work folder and a record are the exception: never edited after their close.
- **Never cited by code.** Code, comments, tests and specs never cite a scaffold: no `.scaffolds/` path, no design or
  source file name, no principle `§` id; a decision is cited by the `_spec` section that holds it. The method's
  `no-design-citations` check fails all three in the code it governs; specs and comments elsewhere are held to it by
  review. Outside `.scaffolds/`, only the root `AGENTS.md`'s table and the governed package's
  `AGENTS.md` point at a scaffold. Another scaffold's designs may cite one.
- **A module moves into the package when a second scaffold needs it.** A scaffold's own principles modules (Seq\*'s
  `queries` and `agent-interface`, from meso\*) stay in its instance until then (the author, 2026-11-23: "keep inside
  for now, once we have more scaffolds we generalize"). Likewise a scaffold's own check moves into `_method/_checks/`
  when a second scaffold needs it, and the rules an instance's fold-backs add stand in its manifest until the package
  takes them.
- **The method learns from every close.** A fold-back's loop changes go into the next loop and into `_method/`'s
  templates, in words with no scaffold in them.
- **Agent context**: never auto-formatted (`.prettierignore`, `biome.json`), machine-invariant (no local paths).
- **`.v6/` is not a scaffold yet**: it moves to `.scaffolds/v6/` in later work, not now (the author, 2026-11-23: "Moves
  later").

## Checks

Each from the repository root, reading the git index (stage first): the package's three (`principles/_checks/`, with
`--package .scaffolds/principles --instance .scaffolds/<name>`), the method's (`_method/_checks/`, with `--instance
.scaffolds/<name>`), and the scaffold's own; each scaffold's `AGENTS.md` lists its commands and arguments.
