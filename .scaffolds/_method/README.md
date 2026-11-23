# The method

How a development scaffold's work is run, written once without any one scaffold in it: templates for the files of a
work folder and a design, and the checks every scaffold runs on itself. It was written from Seq\*'s slices (the
fifth's loop with its fold-back, and the rules the earlier fold-backs added); a scaffold adapts a template where its
work differs, and says so in its loop.

The rules behind it are the principles package's ([`_core`](../principles/_core/README.md), the work-folder files in
[workflow-files](../principles/_core/patterns/workflow-files.md), drift in
[drift-triage](../principles/_core/patterns/drift-triage.md)); what a scaffold is and its rules are in
[../AGENTS.md](../AGENTS.md).

## Templates

| File | For | Copied to |
|---|---|---|
| [loop.md](loop.md) | how the orchestrator runs a work: lanes at once, each step from brief to record, pausing, running alone, gates, checks, stop and ask, finalize | `execution/<work-id>/loop.md` |
| [plan.md](plan.md) | what and why, written once: stages and lanes, gates, the designs and their lanes, measurements as budgets; a re-plan a dated line; no status, no ticks | `execution/<work-id>/plan.md` |
| [checklist.md](checklist.md) | one line per item (lanes, gates, fix lanes, the novice pass, the measurements, FINALIZE), ticked with its commit | `execution/<work-id>/checklist.md` |
| [state.md](state.md) | the current state only, rewritten in place, about a screen: where it stands, in progress, next, blockers, verify | `execution/<work-id>/state.md` |
| [findings.md](findings.md) | append-only: a row per merged lane and gate, the findings by kind | `execution/<work-id>/findings.md` |
| [status.md](status.md) | optional, an instance's routing page: built and open now, edited in place | `status.md` |
| [brief.md](brief.md) | one lane's whole instructions: model, base, numbered items with their proof, where, checks, report | `execution/<work-id>/briefs/<lane>.md` |
| [fixes.md](fixes.md) | the fix phase's queue of small findings that need no design | `execution/<work-id>/fixes.md` |
| [drift-note.md](drift-note.md) | a disagreement between code, spec and design, with options and a recommendation | `execution/<work-id>/drift/<topic>.md` |
| [foldback.md](foldback.md) | what the work taught, the earlier rules checked, the self-check, the author's decision line | `execution/<work-id>/foldback.md` |
| [record.md](record.md) | the permanent record: asked, decided, delivered, learned, left open | `records/<work-id>/<date>-<name>.md` |
| [design.md](design.md) | an agreed design: intent, precedents, decision, defaults with a mock, configuration, consequences, tensions, in the specs | `design/<topic>.md` |

A template's comments say what goes in each section and why; delete them as the section is filled.

The work folder's files follow the package's (§plan-and-checklist, §state-file; the templates in
[workflow-files](../principles/_core/patterns/workflow-files.md)): the plan says what and why and is not edited for
status; the loop reads `checklist.md` and ticks it with each commit; `state.md` is rewritten, never appended (the
closed-work check notes one past 600 words); a lane's and a gate's findings are rows of `findings.md`; Finalize starts
from the checklist's unticked items, and a closed folder has none.

## Checks

In [`_checks/`](_checks/), each run from the repository root with the instance it checks:
`node .scaffolds/_method/_checks/<check>.ts --instance .scaffolds/<name> [flags]`. A check's scaffold-specific part (a
spec map, a baseline, a scope) is a JSON file in the instance's `_checks/`, named in the check's header; a check
whose file is missing says so, or runs with its defaults where it has them. Checks only one scaffold needs stay in
that scaffold's `_checks/`.

| Check | Rejects | Instance file |
|---|---|---|
| [closed-work.ts](_checks/closed-work.ts) | a closed work folder without `foldback.md` and the author's decision line, or with an unticked `checklist.md` line; notes an open folder's `state.md` past 600 words (`--self-check`) | `closed-work.json` (folders closed before the rule) |
| [defaults-designed.ts](_checks/defaults-designed.ts) | a design whose Consequences name something drawn, with no `## Defaults` | `defaults-designed.json` (the words of a drawn thing; exempt designs, only shrinks) |
| [designs-in-specs.ts](_checks/designs-in-specs.ts) | a design with a decision and no `## In the specs` row resolving to a spec section or `design-only: <reason>` | `designs-in-specs.json` (the spec map; pending designs) |
| [no-design-citations.ts](_checks/no-design-citations.ts) | code naming a design or source file of the scaffold, a `.scaffolds/` path or a principle `§<id>` | `no-design-citations.json` (the scaffold's packages, the pointer files) |
| [file-size.ts](_checks/file-size.ts) | a source file past §module-size's bound, or past its baseline | `file-size.json` (the folders), `file-size-baseline.json` |
| [features-planned.ts](_checks/features-planned.ts) | a features row still planned when every lane of its design is merged in the plan | none (`--plan`, the features table at `design/features.md`) |
| [links.ts](_checks/links.ts) | a relative link or a `node .scaffolds/…` command in `.scaffolds/` that does not resolve | none (closed work folders and records reported, not failed) |

They read tracked files (the git index): stage first. They import [`_checks/lib.ts`](_checks/lib.ts), which wraps the
principles package's `lib.ts` (its arguments with `--instance` required, the instance's files); only `links.ts`, which
takes no instance, imports the package's directly.

## A scaffold's life

1. **Start**: ask the author the adopting questions ([adopting](../principles/_core/patterns/adopting.md), step 0)
   with the answers fixed by this layout: the package is `.scaffolds/principles/`, the instance `.scaffolds/<name>/`.
   Create the instance from the package's manifest template: `AGENTS.md` (what it builds, where things are, its
   Checks), `status.md`, `principles/README.md` (the manifest: provenance, folder map, modules, adaptations),
   `_skills/<name>/SKILL.md`, `sources/`, `design/`, `records/`, `execution/`, `_checks/` (the method's instance files,
   and its own checks). Add its pointer to the root `AGENTS.md`'s table and to the governed package's `AGENTS.md`;
   nothing else outside `.scaffolds/` names it.
2. **Design**: drafts in `execution/<work-id>/drafts/`, agreed by the author into `design/<topic>.md` from
   [design.md](design.md), each with its mock when it adds something visible.
3. **Run**: a work folder per slice, from [plan.md](plan.md), [checklist.md](checklist.md), [loop.md](loop.md),
   [state.md](state.md) and [findings.md](findings.md); a brief per lane; the loop runs until its checklist is ticked.
4. **Close**: the loop's Finalize: reconcile, `foldback.md` with the author's decision, the gate, the record, the
   closing line of `state.md`, the next plan. The work folder stays, never moved or deleted unless the author says so.
5. **Re-sync**: when the package changes upstream, see below.
6. **Keep**: when the governed packages stand on their own, the scaffold stays, closed and editable. What the code
   follows has moved into the packages' `_spec/`; the scaffold keeps the reasons, the history and the work, and is
   never cited by code, comments, tests or specs.

## An instance and its package

- **An instance keeps only what is its own.** Its manifest holds its own rules and how it applies a package rule (a
  number, a path, a tool, a tightening), one line each, within the package's rule budget; the evidence behind a line
  is in the records, not the manifest. Why: re-adapted after a re-sync, one instance's reading path fell from about
  34.3k to 26.5k tokens (its manifest from 32k to 9k characters), and the re-sync alone had added about 1k.
- **Updates re-sync from upstream, never fork.** The package folder is replaced by upstream's, edited only as the
  manifest's Provenance says; a change the package needs is proposed upstream, never made in the copy.
- **A rule dropped as "the package's" names what covers it.** Re-adapting an instance after a re-sync, each rule
  dropped because the package states it names the package's id and quotes the clause that covers it, in the
  re-adaptation's report. Where the package states it more weakly ("may", "optional", a default of none), the
  instance keeps one line that tightens it, naming the id. Why: three rules were weakened by being dropped as the
  package's when the package said less (a fold-back at every close, the standing novice pass, a mock's deletion).
