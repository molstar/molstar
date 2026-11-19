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

## In flight when paused: HC (the hover card and the keys panel)

- **Where**: branch `seq/hc-hover-card`, worktree `.claude/worktrees/agent-aea7521d4a8662db8`. The brief is
  `briefs/hc-hover-card.md`; the design is ux-foundations §9, ports.md and viewer-models.
- **Status**: the first round is committed (`9af53e23c`, `4eba79cc8`) and reviewed. Its fix round was sent and was
  running at the pause. The author asked to merge HC and then pause.
- **If the fix round finished**: its commits are on the branch after `4eba79cc8`, starting with a merge of
  `packages-seq` 3aa4db207.
- **If it did not finish**: start a new strong lane on the branch with the list below as its brief.

The fix round's list (the review's findings, as sent):

1. **The port's place** (decided: the author answered "Accept"; the review's top-level `SeqViewerSpec.hoverCard`
   field was offered and not answered). The port stays the `hoverCard` behavior's `port` param.
   - The default spec supplies `hoverCardOverlay` (no hard import in the behavior).
   - An ownership flag like history's and animation's (`disposeWithViewer`).
   - The port is made once, so a show before mount doesn't make a second one at `attach`.
   - The viewer `_spec` says so.
2. **Keys inside the keys panel stay in it**: stop keydown propagation except Escape, or the container's handler
   ignores other targets. Enter on Close closes the panel.
3. **Focus**: the panel takes focus on open, gives it back on close, and handles Esc itself. It is non-modal, and the
   spec says so.
4. **A press clears the rest timer** (`interrupt` and `leave`).
5. **No stale card**: hide it or re-place it at the cursor when the view or the frame changes (zoom, pan, edit, model
   load, a source's labels arriving). The cursor's card is built after the next frame.
6. **`attach`**: an optional `attach(container, viewer)` declared on the behavior type, instead of duck typing.
7. **Esc**: it cancels a running drag before it closes the panel.
8. **`hover.card`**: turned off while shown, it hides the card. Turned on after mount, it works, or the spec says it
   is read at mount.
9. **One `Describer`**: the announcer and the card share it. The card's row elements are reused.
10. **Units and missing labels**:
    - A value shows its unit ("Å" when declared), and Complex's test asserts the exact text.
    - A missing numbering label shows "–".
    - The tautological `said === card.text` check goes, or stays beside the literal checks.
11. **Public API**:
    - Stop exporting `anchorOf`, `hoverCardBehavior`, `keysPanelBehavior`, `mountKeysPanel`, `CardGap`, `valueText`.
    - On `Describer`, only `describe` and `text` are public.
    - `DefaultToolbar`'s `help.keys` goes in the changelog.
12. **Size**: confirm with `pnpm -C packages/seq/core size` that the minimal viewer bundles no card or panel.
13. **Ruler labels** (the author: a ruler in its own numbering needs no label, ux-foundations §9): drop
    `label: 'Position'` in `examples/src/p53/layers.ts` and `trna/trna.ts`. Keep Kabat and IMGT. Regenerate the
    affected records.

**At HC's merge**:
- Run the checks on the merged tree (`loop.md`, Checks), plus folder durations.
- Look at the card in Storybook, light and dark.
- Tick HC in `plan.md` and append `state.md`.
- `.v6/changelogs/seq-hover-card.md`, from the lane's text: the hover card through a port, `hover.card`, the keys
  panel `?` and the toolbar's `help.keys`, `describeKeys`.
- `ports.md`'s table: the hover card port lives on the `hoverCard` behavior's `port` param, not in `spec.ports`. State
  the rule for UI ports a behavior owns.

## Next, in order

1. **RK** (strong): the keyboard cursor's rows, decided "Yes, after HC". Brief `briefs/rk-cursor-rows.md`. Start it
   from `packages-seq` after HC merges.
2. **Stage 7, M3** (orchestrator): the heavy lane (`pnpm test:heavy`, never two heavy runs at once) and the opt-in
   tiers; the browser at dpr 2; PNG export; the bundles.
   - **The minimal viewer fix needs the author's call before a lane**: registries injectable or lazy by default, lazy
     imports, `sideEffects: false`. It changes defaults and public API.
   - Then the GATE: the author confirms or moves the targets.
3. **The use gate** (deferred): the author uses the views. At finalize, only a look at what changed since.
4. **FINALIZE**: reconcile designs, `_spec`s, status, features, edge cases, `TODO.md` and `drift/`. The fold-back is
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
