# Plan: the sixth slice (draft)

**Draft, 2026-11-23**: not agreed and not started. It starts after the fifth slice closes (its stage 5 and finalize),
and runs alone: it moves every path the loop, the briefs and the checks use, so no lane may be running.

Work id `sixth-slice`. **Scaffolds**: `.seq/` reconciled into a common place for this and future development scaffolds,
instead of deleted.

- The author, 2026-11-23: "originally i planned to remove it, but i think we should reconcile it instead and put it in
  some common place, where next similar "scaffolds" will live".
- On the proposal of a common folder: "let's do the scaffolds as 6th slice. do a draft plan there for this. also, "oz"
  was a temporary working name. maybe we don't even need to use it and just do "principles" etc".

## What `.seq/` holds and where each part goes

| Part | Now | Kind | Goes to |
|---|---|---|---|
| The principles package | `.seq/oz/` (copied from the oz exploration repository at `7d03f60`) | generic, copied | `.scaffolds/principles/`, renamed: no "oz" in names or text, its provenance in one line |
| Its skills | `.seq/oz/_skills/{principles,principles-router,principles-self-check}` | generic | with the package |
| Its checks | `.seq/oz/_checks/` (ids, paths, merge gate) | generic | with the package |
| The execution method | the slices' `loop.md`, briefs, `state.md`, `plan.md`, `fixes.md`, `foldback.md`, the record's shape | generic in shape, Seq\*'s in text | `.scaffolds/_method/` templates, written from the fifth slice's files with Seq\* taken out |
| The instance's checks | `.seq/_checks/` (file size, EC ids, features-planned, defaults-designed, closed-work) | half generic | generic ones parameterised and moved to the method; Seq\*'s baselines stay in its instance |
| The modules from meso\* | `.seq/principles/{queries,agent-interface}/` | generic rules from another instance | a question (below): into the package, or kept in the instance |
| Seq\*'s manifest and skill | `.seq/principles/README.md`, `.seq/_skills/seq/` | Seq\*'s | `.scaffolds/seq/` |
| Designs, sources, edge cases | `.seq/design/`, `.seq/sources/` | Seq\*'s knowledge, with its reasons and history | `.scaffolds/seq/`, read-only, after the reconcile below |
| Records | `.seq/records/` | Seq\*'s history, small | `.scaffolds/seq/records/` |
| Execution folders | `.seq/execution/` (2.3 MB) | Seq\*'s history, large | a question (below): kept in the tree or left to git |

## Layout (proposed)

```
.scaffolds/
  AGENTS.md        what a scaffold is, how one is started, run, closed and archived
  principles/      the principles package (was .seq/oz)
  _method/         templates: loop, brief, state, plan, fixes, record; generic checks
  seq/             Seq*'s instance: manifest, skill, design/, sources/, records/ (execution/ if kept)
```

A next scaffold starts as `.scaffolds/<name>/` from `_method/`, its manifest listing the package's modules and its own.

## Stages

- [ ] **Stage 0: decisions** (the questions below), and this plan agreed.
- [ ] **Stage 1: reconcile the specs**. Every design's decisions that the code follows are in a package's `_spec/`
      (core, viewer, ext-rondo, examples where it has one): a pass per design, its gaps written into the spec. Done
      when a check lists no design decision without a spec section citing it, or a reason it stays design-only.
- [ ] **Stage 2: rename oz to principles**, in place first (`.seq/oz` → `.seq/principles-package` or straight to the
      target in stage 3): folder names, text ("the oz principles" → "the principles"), skill descriptions, the
      manifest's folder map, the checks' `--package` argument. One provenance line keeps the upstream repository and
      commit. A check that no "oz" is left outside that line.
- [ ] **Stage 3: move** with `git mv` (history kept): the layout above. Rewrite every path that names `.seq/`: the
      instance's links, the checks, the skill, `packages/seq/AGENTS.md`'s one pointer, `biome.json`,
      `.prettierignore`, `.gitignore`. Confirm nothing in code cites `.seq/` (the rule; `cli/seqx-render` and
      `cli/structure-info` match `.seq` in their text: check that it is a file extension, not a link).
- [ ] **Stage 4: the method**: templates in `_method/` from the fifth slice's files (loop with its fold-back, brief,
      state, plan, fixes, record), the generic checks taking the instance as an argument; `.scaffolds/AGENTS.md`.
- [ ] **Stage 5: verify**: every check run from the new place on the Seq\* instance; a link check over `.scaffolds/`
      (no dangling relative links); the root `AGENTS.md`'s table gains the `.scaffolds/` row; `.seq/` gone.
- [ ] **Docs** (the author, 2026-11-23: "we should write "Docs" folder in storybook in slice 6 too, that explains
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

## Questions for the author (stage 0)

1. **The name**: `.scaffolds/` (proposed), or another.
2. **The execution folders**: kept in the tree under `seq/execution/` (browsable, 2.3 MB), or left to git with the
   records as the index (each record names its slice's last commit).
3. **The meso\* modules** (`queries`, `agent-interface`): promoted into the package (they are generic rules), or kept
   as Seq\*'s own modules.
4. **`.v6/`**: stays where it is (a release line's designs, not a scaffold), or becomes `.scaffolds/v6/` later.
5. **`.seq/AGENTS.md`'s rule "temporary, deleted"** is replaced by "archived, read-only, never cited by code": the
   instance stays citable from other scaffolds' designs, never from code or specs.
