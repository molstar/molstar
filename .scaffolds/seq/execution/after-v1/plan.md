# Plan: after version 1

**Draft, 2026-11-24**: written at the close-out's close (its loop, Finalize step 5), not agreed and not started. The
work id `after-v1` is a placeholder; the author names the milestone at its Gate 0.

Version 1 of Seq\* (next) is closed ([the record](../../records/close-out/2026-11-24-close-out.md)): the language at
`0.11.0`, every planned row built, deferred with its trigger, or cut. What is left is carried work, the open questions
the close-out found, and the step from `packages-seq` to the repository's main line.

## Where this is heading

| Work | Milestone | In one line |
|---|---|---|
| first to sixth slice | the language, the viewer, genomes, the ring, the scaffolds | closed, each with its record |
| close-out | version 1 | closed 2026-11-24 |
| after-v1 | carried work and the merge to the main line | this draft |

## Carried from the close-out

- **The 200-copy split and merge over their edit bound** (the author: "Carry it"; `packages/seq/TODO.md`): the split's
  first frame 21.1 ms and its task 28.7 ms, the merge's task 17.5 ms. Profile the remaining microtask tails first; then
  a design question: lays shared between identical copies.
- **The novice pass's looks and small gaps** (`packages/seq/TODO.md`, its close-out line): a partner step on the ring,
  search's input hints, the Ligands story's legend, a group's fold in history, the looks C6 to C17, B5 to B13, A1 to
  A12.
- **The square's whole genome after a ring cut** and **a header's pressed state** (`packages/seq/TODO.md`).
- **Seven edge cases not met** (EC-34, 52, 98, 104, 105, 158, 168 in part; `design/edge-cases.md`, Outcomes of the
  close-out): each built, or its row's status says why not.
- **The lesser measurement overs** (`packages/seq/TODO.md`): the stories' settled frames, the inferred GC spikes.
- **Outside Seq\* (next)**: the old Seq\* viewer's label churn over its allocation baseline (`packages/sequences/viewer`,
  its owner's).

## Gate 0 (questions for the author)

1. **The milestone**: the merge of `packages-seq` to the main line first (the PR, its review, what a reviewer is
   pointed at), or a carried-work slice first? Recommended: the merge first, the carried work after it on the main
   line.
2. **The split patch's design**: shared lays between identical copies now, or deferred until a dataset of copies
   needs it? Recommended: profile first, decide on the numbers.
3. **The deferred rows** (keyed tables, keep this region, map kinds, A3M and Clustal, bigWig's own levels, the glyph
   atlas): any trigger met? Recommended: none yet.

## How this work runs

The close-out's loop with its fold-back applied (B1, B2, B4, B5: `../close-out/loop.md`), three lanes at once, every
lane's look in a visible pane before its merge.

## Done when

Written at Gate 0 with the milestone.
