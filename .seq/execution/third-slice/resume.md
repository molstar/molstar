# Resume: the third slice, paused 2026-11-19

This file is for an orchestrator picking the slice up in a new session. Read `loop.md` (the rules), then `state.md`
(its tail), `plan.md` (the unticked items), `gaps.md`, `followups.md` and `drift/`. Restart with
`/loop run .seq/execution/third-slice/loop.md`.

## Where it stands

- **Branch** `packages-seq` (main checkout `/Users/davinci001/projects/molstar/molstar-v6`). Every lane merges here
  `--no-ff`; nothing is pushed.
- **Merged this slice**:
  - files and records: FS-A1, FS-A2;
  - looks and fixes: LG, LS, DX, HB, CR, IZ, HL, PR, SW, AY, OP;
  - rows and queries: RC, LB, QS;
  - datasets: DP, DT, DG, DC, DA;
  - language and core: NB (numbering kind), CY (source-backed axes, `openModel`), TW (track padding, squircle
    wells and emphasis);
  - tests and the novice pass: QT2 (quick-lane time), NP (the novice pass's fixes).
- **Gates passed**: stages 1–5. Stage 6: AY, OP, SW and the novice pass done (`novice-findings.md`); NP merged.
- **Last checks on `packages-seq`** (after NP): quick 2585, medium 105, perf 127; tsc, lint, format and workspace
  green. Sizes (after DA): reading 211.2, headless 676.7, viewer default 1063.6, minimal 999.6 kB (all over the
  guides: stage 7).
- **Quick lane time**: at QT2's merge, core 1.98, viewer 1.96, examples 1.08, storybook 1.96 s. Lanes since then
  measured 2.4–3.5 s under load 6–13, which is not comparable. Re-measure on a quiet machine (`uptime` under 3) before
  the next gate.

## HC: merged before the pause

HC merged (`16112a9b7`) after its review round; recorded (`.v6/changelogs/seq-hover-card.md`, ports.md's rule for a UI
port a behavior owns). On the merged tree: quick 2606, medium 105, perf 127, workspace green; sizes viewer default
1085.5, minimal 1010.9 kB (the minimal bundles no card or panel, an esbuild metafile checked; the size script's
`--files` crashes after its first entry). Quick lane at load 4.2: core 2.33, viewer 2.53, examples 1.29, storybook
2.46 s: over the 2 s guide; re-measure quiet (load under 3) first, then a QT3 lane if still over. HC's leftovers:
rulers naming their chain (Complex's, Antibody against's) kept labelled; an empty track label still drawn as empty
text (core `tracks.ts`); a cursor card replaced when a zoom re-hits the hover.

## Next, in order

1. **The quick lane's time**: re-measure quiet; QT3 (fast) if over the 2 s guide.
2. **RK** (strong): the keyboard cursor's rows, decided "Yes, after HC". Brief `briefs/rk-cursor-rows.md`. Start it
   from `packages-seq` after HC merges.
3. **Stage 7, M3** (orchestrator): the heavy lane (`pnpm test:heavy`, never two heavy runs at once) and the opt-in
   tiers; the browser at dpr 2; PNG export; the bundles.
   - **The minimal viewer fix needs the author's call before a lane**: registries injectable or lazy by default, lazy
     imports, `sideEffects: false`. It changes defaults and public API.
   - Then the GATE: the author confirms or moves the targets.
4. **The use gate** (deferred): the author uses the views. At finalize, only a look at what changed since.
5. **FINALIZE**: reconcile designs, `_spec`s, status, features, edge cases, `TODO.md` and `drift/`. The fold-back is
   its own step, decided by the author. Write the record and the fourth slice's plan.

## Open questions for the author

- `gaps.md`, Open: the lists-free palette module, colour-blind-safe substitutes, dim dark letter cells, and the rest
  listed there.
- **`itWithin` and the medium lane**: should `itWithin` gain an option that puts its unbounded half in medium? QT2
  wrote the plasmid pan by hand. `@molstar/test-harness` is shared.
- **The minimal viewer** (stage 7, above).
- **ANARCI**: record its numbering of the Antibody design by hand, to check the design's two rules (DA's followup).

## Decisions taken this session

Each is also in `state.md`, with the author's words:
- The hover card and keys panel this slice (HC), through a port.
- Every emphasis is a squircle.
- Track padding and squircle wells (TW).
- CY's lengths through resolve, then `openModel`.
- The numbering kind in 0.4.0 (NB).
- Alignment-column pairs, ranges and maps go to 0.5.0, in slice 4.
- A ruler in its own numbering is unlabelled.
- The cursor's rows after HC (RK).

## Housekeeping

- **`.claude/launch.json`** has a local, uncommitted `storybook-lane` entry (port 6008) pointing at the main
  checkout. Point it at a lane's worktree to look at that lane's stories. Revert it at finalize; never commit it.
- **A Storybook server may still be running on 6008**, from the old session's browser pane. Stop it, or reuse it.
- **Old lane worktrees** remain under `.claude/worktrees/` (their branches are merged).
  - Clean them at finalize with `git worktree remove`, after checking each branch is merged into `packages-seq`.
  - `seq/da-antibody` (worktree `agent-ac15be1344bf444d0`) was never used. The work merged from `seq/da-antibody-b`.
- **The record script**: `node --experimental-strip-types packages/seq/examples/scripts/record.ts <dataset>`, after
  `npx tsc -b packages/seq/examples`. It fetches the network: the author's approval per run.
