# Gaps: close-out

Decisions a gate asks for (the loop, Gates).

- **Decided** (the author, 2026-11-24: "you can increase the budget if needed"): raised to 1,095 kB with its cause. **The minimal viewer's budget**: 1,089.6 kB against the 1,090 the author accepted for version 1 (2026-11-23), after
  CL1 and CL4; 0.4 kB left, and wave 2 and 3 still draw. Recommended: lazy-load the region-read and blocks code (about
  60 kB, the stage 0 option B) as a lane before CL6, rather than raise the budget. Asked at GATE 1.
- **Decided** (the author, 2026-11-24: "ok with the right click"): kept, written into the design and spec. **A right-click on a ring's ribbon**: the ribbon note (resolved, deleted). Asked at GATE 1.
- **The keys for the next feature**: decided, `.` and `,` (the author, 2026-11-24); to CL3.

## Decided overnight

Each made as recommended while the author was away (2026-11-24 night), for the morning's reconcile.
- **The zoom hatch over a track that draws nothing** ([drift/zoom-hatch-empty.md](drift/zoom-hatch-empty.md)): the hatch where data exists but this level does not draw it; the empty-well row where the region has none. The author to confirm.
- **Search (CL5) and track groups (CL6) have no design**: their lanes build a core and need one first. Overnight, two
  drafts are written (`drafts/search.md`, `drafts/track-groups.md`, no code) for the author's agreement in the morning;
  CL5 and CL6 start after it. Wave 3 runs CL9 first (after CL3, which shares its cursor code), then wave 4's lanes
  with designs (CL7, CL8, CL13) and FX25, before CL5 and CL6. Recommended: agree the drafts, or amend them.
- **Track groups need a language change** (`drafts/track-groups.md`: a `group` node, 0.10.0 → 0.11.0, no migration):
  a language change is never decided overnight; CL6 waits for the author's word on the draft.
- **Search's draft** (`drafts/search.md`) adds public API (a matcher registry, `SearchModel`, `search.*` commands,
  two highlight channels, Mod+G and F3): within what the author agreed for search, but an API is never decided
  overnight; CL5 waits for the author's word.
- **FX25 stopped on a permission refusal** (2026-11-24 night; the loop's stop condition, never worked around): the
  `cc-plugin-privacy` plugin refused two of the lane's Bash commands with "it, or a script it names, holds or reads the
  account's email address beside what looks like a client's identity, or puts the address in a URL". Neither command
  held an email; both named Seq\* source files, whose headers carry the author's `@author` line (likely the match). The
  refused commands: (1) a `python3 -` script editing `core/src/render/lanes.ts` and `ranges.ts` (the lane then made the
  same edits with the Edit tool: past the stop rule, reported by the lane itself); (2) a `cp`/`git checkout`/vitest run
  showing item 1's test failing on its parent. FX25's worktree (`.claude/worktrees/agent-a98a4bf365b8a5323`, `seq/fx25`)
  holds items 1 and 2 uncommitted (TP53-201's lane: sticky lanes repacked at settle for a marked item, a 146-view
  canvas-versus-SVG lane sweep; a variant's card titled by its ID); items 3–6 not started. **For the author**: allow or
  adjust the plugin for such commands, or say how FX25 goes on; the orchestrator does not run the refused commands.
- **FX25's question**: an allele's card on the Sites track keeps its spec'd title (`pos:ref>alt`, "never its record's
  ID"), while a genotype cell's card now titles by the record's ID. Recommended: keep both (the spec's rule for alleles).
- **`rowsValue` over a lazy alignment needs a data design** (CL13 item 4, stopped): a lazy source serves level 1 only
  and no row summary, so one value a row needs a sampled one-row read per device row, or a source-side row summary
  (tile budgets, `layerRefused`, the non-gap share only). A new data design is not decided overnight. Recommended:
  carry it past version 1 (the TODO's Navigators line already holds it), the navigator still saying "a lazy alignment's
  rows are left out".
- **A drag's codon snap over a codon split by an intron** (CL9's review): decided overnight as recommended: the split
  codon is selected as its pieces, never the intron's nucleotides (the design draws a split codon in each piece); the
  code's "hull" changed to that. The author to confirm.
- **Values marks (CL8) need a language change** (`drafts/values-marks.md`): points (with the min-to-max column) and
  the log scale are built already (the fifth slice's VP, EC-215); what row 130 still lacks is reference lines: a
  `reference` param on the values layer (up to 4 lines, `at`, `label`, `color`, `style`), a minor bump. Never decided
  overnight. The two language changes (track groups' `group` node, values' `reference` param) would share one bump if
  both are agreed: one language lane, in order. Recommended: agree both drafts, one language lane bumping to 0.11.0.

- **A 0.x version's minor read as a newer minor** (CL9, the version-aware validate): a 0.11.0 spec's unknown kind is a
  warning to a 0.10.0 reader, as the brief decided, though 0.x versions usually break on a minor. Recommended: keep it
  (the language is ours, and an older reader drawing what it knows is the point); revisit at 1.0.
