# Plan: the sixth slice (draft)

**Draft, 2026-11-23**: not agreed and not started. It starts after the fifth slice closes (its stage 5 and finalize),
and runs alone: it moves every path the loop, the briefs and the checks use, so no lane may be running.

Work id `sixth-slice`. **Scaffolds**: `.scaffolds/seq/` reconciled into a common place for this and future development scaffolds,
instead of deleted.

- The author, 2026-11-23: "originally i planned to remove it, but i think we should reconcile it instead and put it in
  some common place, where next similar "scaffolds" will live".
- On the proposal of a common folder: "let's do the scaffolds as 6th slice. do a draft plan there for this. also, "oz"
  was a temporary working name. maybe we don't even need to use it and just do "principles" etc".

## What `.scaffolds/seq/` holds and where each part goes

| Part | Now | Kind | Goes to |
|---|---|---|---|
| The principles package | `.scaffolds/principles/` (copied from the oz exploration repository at `7d03f60`) | generic, copied | `.scaffolds/principles/`, renamed: no "oz" in names or text, its provenance in one line |
| Its skills | `.scaffolds/principles/_skills/{principles,principles-router,principles-self-check}` | generic | with the package |
| Its checks | `.scaffolds/principles/_checks/` (ids, paths, merge gate) | generic | with the package |
| The execution method | the slices' `loop.md`, briefs, `state.md`, `plan.md`, `fixes.md`, `foldback.md`, the record's shape | generic in shape, Seq\*'s in text | `.scaffolds/_method/` templates, written from the fifth slice's files with Seq\* taken out |
| The instance's checks | `.scaffolds/seq/_checks/` (file size, EC ids, features-planned, defaults-designed, closed-work) | half generic | generic ones parameterised and moved to the method; Seq\*'s baselines stay in its instance |
| The modules from meso\* | `.scaffolds/seq/principles/{queries,agent-interface}/` | generic rules from another instance | kept in Seq\*'s instance for now (decided, 1 below) |
| Seq\*'s manifest and skill | `.scaffolds/seq/principles/README.md`, `.scaffolds/seq/_skills/seq/` | Seq\*'s | `.scaffolds/seq/` |
| Designs, sources, edge cases | `.scaffolds/seq/design/`, `.scaffolds/seq/sources/` | Seq\*'s knowledge, with its reasons and history | `.scaffolds/seq/`, read-only, after the reconcile below |
| Records | `.scaffolds/seq/records/` | Seq\*'s history, small | `.scaffolds/seq/records/` |
| Execution folders | `.scaffolds/seq/execution/` (2.3 MB) | Seq\*'s history, large | `.scaffolds/seq/execution/`, kept in the tree (decided) |

## Layout (proposed)

```
.scaffolds/
  AGENTS.md        what a scaffold is, how one is started, run, closed and archived
  principles/      the principles package (was .scaffolds/principles)
  _method/         templates: loop, brief, state, plan, fixes, record; generic checks
  seq/             Seq*'s instance: manifest, skill, principles/ (meso*'s modules), design/, sources/, records/, execution/
```

A next scaffold starts as `.scaffolds/<name>/` from `_method/`, its manifest listing the package's modules and its own.

## Stages

- [x] **Stage 0: decisions** (2026-11-23, below). The plan itself is agreed when the slice starts.
- [x] **Stage 1: reconcile the specs** (done 2026-11-23: S1A `936242251`, S1B `81a82d453`, S1C `3001f67d2`, S1D; 46 designs placed, `PENDING` empty; the gate green: quick 4,887, lint, format, every scaffold check; 15 drift notes to the author) (started 2026-11-23: lanes S1A–S1D, `briefs/s1-spec-pass.md`; the check is
      `.scaffolds/seq/_checks/designs-in-specs.ts`, done when its `PENDING` is empty). Every design's decisions that the code follows are in a package's `_spec/`
      (core, viewer, ext-rondo, examples where it has one): a pass per design, its gaps written into the spec. Done
      when a check lists no design decision without a spec section citing it, or a reason it stays design-only.
- [x] **Stage 1b: comments cite specs, not designs** (done 2026-11-23, S1E `5b09abdaa`; `.scaffolds/seq/_checks/no-design-citations.ts`) (found 2026-11-23 during stage 1: about 160 comments and test
      headers in `packages/seq` cite a design by file name, e.g. "molecules.md, Links"; a design is the scaffold's). A
      lane (Sonnet; mechanical) rewrites each to the spec section the design's `## In the specs` table names for that
      decision (`core _spec, Links`), or drops the citation where the row is design-only; a package's own `_spec/*.md`
      (e.g. core's `util/_spec/assets.md`) stays. Done when a grep for the designs' file names in `packages/` finds none
      (a check, kept with the method's checks). Also: comments citing a spec section by an old name (`data/_test/bounds.test.ts`, "Data, Bounds", renamed "An answer's bounds" by S1C), and the core spec's own citations of design files (about line 587). Runs after stage 1, before stage 3.
- [x] **Stage 2: rename oz to principles** (done 2026-11-23, `0be663880`: the live text; closed slices and records keep their history), in place first (`.scaffolds/principles` → `.scaffolds/seq/principles-package` or straight to the
      target in stage 3): folder names, text ("the oz principles" → "the principles"), skill descriptions, the
      manifest's folder map, the checks' `--package` argument. One provenance line keeps the upstream repository and
      commit. A check that no "oz" is left outside that line.
- [x] **Stage 3: move** (done 2026-11-23: `88f55a78c` the moves, the next commit the rewrites; every check green from its new place; one link the move broke, fixed; 87 links that dangled before the move, all in closed drafts, left as history) with `git mv` (history kept): the layout above. Rewrite every path that names `.scaffolds/seq/`: the
      instance's links, the checks, the skill, `packages/seq/AGENTS.md`'s one pointer, `biome.json`,
      `.prettierignore`, `.gitignore`; and the five tests in `packages/seq/core` whose header comment cites `.scaffolds/seq/sources/seqstar/edge-cases.md` (found 2026-11-23: the citation goes, the `EC-<n>` id stays). Confirm nothing in code cites `.scaffolds/seq/` (the rule; `cli/seqx-render` and
      `cli/structure-info` match `.seq` in their text: check that it is a file extension, not a link).
- [ ] **Stage 4: the method**: templates in `_method/` from the fifth slice's files (loop with its fold-back, brief,
      state, plan, fixes, record), the generic checks taking the instance as an argument; `.scaffolds/AGENTS.md`.
- [ ] **Stage 5: verify**: every check run from the new place on the Seq\* instance; a link check over `.scaffolds/`
      (no dangling relative links); the root `AGENTS.md`'s table gains the `.scaffolds/` row; `.scaffolds/seq/` gone.
- [x] **Docs** (merged 2026-11-23, `c9d120931`; taken early as lane DX, 2026-11-23, `briefs/dx-storybook-docs.md`: written from the current `_spec`s;
      stage 1's reconciliation checks the docs too) (the author, 2026-11-23: "we should write "Docs" folder in storybook in slice 6 too, that explains
      this very well right of the bat, so people dont expect opiniated viewer like jalview out of this. and in general
      the docs should go high level around how this works etc. similar to what the docs for the previous version do"):
      a `Docs/` folder in the Seq\* (next) storybook, high level, in the shape of the previous version's
      (`packages/sequences/storybook/src/_stories/docs/`: introduction, concepts, principles, format, features,
      extending, query language, streaming, colour schemes), written new for `packages/seq`. It opens with what Seq\* is
      and is not: a library that makes the pictures, the engine and not the game; composition and navigation are the
      application's (the ring and its square as the worked example of an application built on the library). Then how
      it works: the language and its levels, data and the store, layers and layouts, navigators, hover and commands,
      extensions, the CLI and headless use; each page with live examples and links into the stories. Written from the
      reconciled `_spec`s (stage 1), so the docs and the specs say the same thing.
- **Deferred: Prism on the new Seq\*, then removing `packages/sequences` and renaming `seqx` to `seq`** (the author,
  2026-11-23, on the draft's second frame tick: "let's not do slice 7 yet at all then for now; we will do that part
  later, because i was planning to rewrite whole prism anyway. so let's just remove slice 7 for now and we will do it as
  a separate piece of work."). Not this slice. The draft stays as a reference for that work:
  [drafts/prism-on-new-seq.md](drafts/prism-on-new-seq.md) (a Mol\* adapter of the seq data port, the mapping, the
  gaps G1–G7, the questions).
- [ ] **Record** in `.scaffolds/seq/records/sixth-slice.md`.

## Decided at stage 0 (the author, 2026-11-23)

1. **The name**: `.scaffolds/` (as recommended).
2. **The execution folders**: kept in the tree, under `.scaffolds/seq/execution/` (as recommended).
3. **The meso\* modules** (`queries`, `agent-interface`): "keep inside for now, once we have more scaffolds we
   generalize. maybe record that in readme in the .scaffolds folder?". They stay in Seq\*'s instance, and
   `.scaffolds/AGENTS.md` (stage 4) records the rule: a module moves into the package when a second scaffold needs it.
4. **`.v6/`**: "Moves later": it becomes `.scaffolds/v6/` in later work, not in this slice.
5. **The rule "temporary, deleted"** becomes: a closed scaffold is kept and never cited by code, comments, tests or
   specs; its files may still be edited when needed ("i think we can still edit if needed (and prolly will, especialyl
   some of the designs)"). Not read-only. Other scaffolds' designs may cite it.
