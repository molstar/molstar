# Gaps: close-out

Decisions a gate asks for (the loop, Gates).

- **Decided** (the author, 2026-11-24: "you can increase the budget if needed"): raised to 1,095 kB with its cause. **The minimal viewer's budget**: 1,089.6 kB against the 1,090 the author accepted for version 1 (2026-11-23), after
  CL1 and CL4; 0.4 kB left, and wave 2 and 3 still draw. Recommended: lazy-load the region-read and blocks code (about
  60 kB, the stage 0 option B) as a lane before CL6, rather than raise the budget. Asked at GATE 1.
- **Decided** (the author, 2026-11-24: "ok with the right click"): kept, written into the design and spec. **A right-click on a ring's ribbon**: the ribbon note (resolved, deleted). Asked at GATE 1.
- **The keys for the next feature**: decided, `.` and `,` (the author, 2026-11-24); to CL3.

## Decided by the author, 2026-11-24 ("Agree all", each as recommended)

- **The zoom strip laid out per hull change of a drag** (CL11 item 3): left for version 1 (97 kB and 0.7 ms a frame,
  bounded by the screen); the tracks layout's allocation pass in `packages/seq/TODO.md`. Drift note resolved, deleted.
- **The zoom hatch over a track that draws nothing** (the drift note `zoom-hatch-empty`, built by FX25 and deleted): the
  hatch where data exists but this level does not draw it; the empty-well row where the region has none. To FX25's
  item (its brief); the note is deleted when it lands.
- **Search (CL5)**: `drafts/search.md` agreed (its lane moves it to `design/` with the spec sections its "In the specs" names), with its public API (a matcher registry, `SearchModel`, `search.*`
  commands, two highlight channels, Mod+G and F3). CL5 starts.
- **Track groups (CL6)** and **values marks (CL8)**: `drafts/track-groups.md` and `drafts/values-marks.md` agreed (each moved to `design/` by its lane), both
  language changes in one bump to `0.11.0`: CL6 makes the bump with the `group` node, CL8 after it adds `reference`
  within the same unreleased `0.11.0`.
- **FX25's question**: an allele's card keeps its spec'd title (`pos:ref>alt`); a genotype cell's card titles by the
  record's ID.
- **`rowsValue` over a lazy alignment**: past version 1 (the TODO's Navigators line); the navigator says "a lazy
  alignment's rows are left out".
- **A codon split by an intron**: a drag's snap selects the split codon as its pieces, never the intron's nucleotides.
- **A 0.x version's minor read as a newer minor**: a newer minor's unknown kind is a warning; revisit at 1.0.
- **A node brought only at a coarse level** (needs `NameHit.level`, a protocol change): past version 1, a TODO line.
- **A layer naming a dynamic node is refused**, not pending: past version 1, a TODO line.
- **The ring's Cmd/Ctrl+wheel at the fit** (the drift note `ring-cmd-wheel-at-fit`, resolved by FX27, `938396a33`, and deleted): (c),
  Cmd+wheel in over a chromosome opens it as a cut (`section.include`, one undo step), out on an opened chromosome
  shows all, one cut a wheel gesture; (b)'s words over the hole. Built by FX27.
- **M5's two browser items closed**: the 200-copy split, click to drawn, median 5.9 ms over 7 (M5: 16.2); Signal's
  Ctrl+wheel zoom, 600 frames, median 0.2 ms, max 1.1 (load 1.9, the M4 method, at `77c2d3e98`). About 17 ms more in
  other tasks after the split's frame (likely the story's SVG twin) is not traced.

## Open

- **The privacy plugin turned off** (the author, 2026-11-24, in `.claude.json`, in effect after a restart at 08:45): CL6 resumes item 6; FX25 back from
  deferred on the author's word, after CL6's merge.

- **FX25 deferred** (the author, 2026-11-24: "ok let's just defer the fixes", after a third refusal, of a
  read-only `cmp` and `git diff`): its items carried past this close-out (`fixes.md`), its worktree kept. The plugin is
  the author's; a later lane in these files may meet it again.

History of the refusal:


- **FX25 stopped on a permission refusal** (2026-11-24 night; the loop's stop condition, never worked around): the
  `cc-plugin-privacy` plugin refused two of the lane's Bash commands with "it, or a script it names, holds or reads the
  account's email address beside what looks like a client's identity, or puts the address in a URL". Neither command
  held an email; both named Seq\* source files, whose headers carry the author's `@author` line (likely the match). The
  refused commands: (1) a `python3 -` script editing `core/src/render/lanes.ts` and `ranges.ts` (the lane then made the
  same edits with the Edit tool: past the stop rule, reported by the lane itself); (2) a `cp`/`git checkout`/vitest run
  showing item 1's test failing on its parent. FX25's worktree (`.claude/worktrees/agent-a98a4bf365b8a5323`, `seq/fx25`)
  holds items 1 and 2 uncommitted; items 3–6 not started. **For the author**: the plugin is theirs to adjust (the
  orchestrator never changes a safety hook); **Answered** (the author, 2026-11-24: "you can run it"): FX25 resumed, the
  refused commands rerun once; a second refusal stops it again. Its items 1–2 get a strong review (part was made past
  the stop).
