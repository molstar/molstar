# Plan: version 1's close-out

**Draft, 2026-11-23**: written at the fifth slice's close (its loop, Finalize step 5), not agreed and not started. It
runs after the sixth slice (scaffolds), which moves every path the loop uses. The work id is a placeholder; the sixth
slice's layout names its folder.

**Goal**: version 1 closed. Its planned rows that no milestone took are built or cut, the fifth slice's carried items
are done, and the author's use of the genome views and the ring has passed.

## Carried from the fifth slice

- **The use gate** (the author, 2026-11-23: "Carry to next plan"): the author tries the genome views and the ring.
  Under the fold-back's B10, a gate the author plays: each view's story is checked against its agreed picture.
- **The novice passes' 14 design questions** (`../fifth-slice/novice-findings.md`, items 5–13, 21, 25, 26, 28, 33,
  35). Each needs a design answer before a lane:
  - a legend entry hiding its layer;
  - the square's hint and its whole-genome control;
  - synteny's legend;
  - the stem menu's partners and its right-click route;
  - a key for the next feature;
  - a canonical transcript;
  - the VCF's sample names;
  - the zoom-in hint;
  - the stem card's first row;
  - a pair card's second end;
  - Hi-C cards in bases;
  - the links-pending state.

  Under B10 each comes with a mock.
- **The minimal viewer's growth** (+150.7 kB over the fifth slice, `../fifth-slice/measurements.md`): accept it as the
  budget, or plan a trim (lazy chunks for the genome axis and the navigator over blocks first).
- **The drag rule** (`../fifth-slice/fixes.md`): a drag's work bounded by the screen, not by the selection. Its
  instances are the PAE drag's card and the 3,000-contig drag. It is a strong lane with a perf-tier bound.
- **The stage 1 novice decisions F3, F9, F11, F14** (`../fifth-slice/fixes.md`).
- **The fold-back's code changes**: `features-planned` extended to a brief's numbered items (B3); the heavy lane's
  tight bounds isolated in `vitest.config.ts` (B7: the interval lookup, the pairs pyramid's slices, the field theme's
  1.1× ratio, each failing in the full lane and passing alone).
- **Feature rows out of `deferred`** (the author's call): two assemblies on one axis, layers across sections, variants
  by sample, gene hierarchies, the name index. All five are built (the close's reconcile has the evidence).
- **M5's leftovers** (`packages/seq/TODO.md`): the 200-copy split patch over 16 ms, one Signal zoom frame of 25.7 ms,
  Cmd+wheel on the ring at the fit.

- **The sixth slice's drift, decided as recommended** (the author, 2026-11-23; the notes in
  `../sixth-slice/drift/`), the code items:
  - the theme picker's BLOSUM62 text (`theme/provenance.ts`) says what the scheme does: it compares with the
    reference row;
  - a folder-graph test in `packages/seq/core` (each folder an entry, no cycles);
  - a lanes band's history step keyed by item, not `bandTop`;
  - `rowsValue` drawn by a navigator over a lazy alignment (the TODO's Navigators);
  - the info issue `seq.render.pairs.filter` (a data filter not applied to bins);
  - the version-aware validate (an unknown kind of a newer minor a warning), before 1.0.

## Version 1's rows no milestone took

To be listed from `design/features.md` (every `planned` row) and `design/version-1.md` at this plan's stage 0, each
built, cut, or deferred past version 1 by the author's word. Known now:
- the Manhattan family and `split`;
- reference lines;
- keyed joins (EC-224);
- the `pairwise` and `alignment` map kinds;
- the query forms `search`, `onto` and `group`;
- drags snapped to codons;
- track groups and summaries;
- 2bit, BigBed and PAF (designed only).

## Lanes (written at stage 0, the author, 2026-11-23: "Write them in")

The answers the lanes build are [drafts/novice-questions.md](drafts/novice-questions.md) (A to G) and
[drafts/version-1-rows.md](drafts/version-1-rows.md). Each answer enters its design, its page and its spec with its lane.

| Lane | Items | Model |
|---|---|---|
| CL1 Legend and layers | A (5, 13, look 24); E7 synteny's legend; legends state their domain (113) | strong |
| CL2 Ring | B (8, 26, 28, 29); D21; D25; G35 on the ring | strong |
| CL3 Square and keys | C (6); D33; E9; G12; G35 on the square | Sonnet |
| CL4 Formats | F10 a canonical transcript; F11 the VCF's sample names; look 34 | Sonnet |
| CL5 Search | the matchers (exact, degenerate, regex), bounded; the viewer's search state and commands, a highlight channel; a matcher registry; a story's search box (182, 180, 181, 207) | strong |
| CL6 Track groups | groups, summaries and collapse, EC-103 (159, 70, 141, 185); the Ranges' collapse (128) | strong |
| CL7 Streamed patches | a bookmark or ref naming a dynamic node, waiting while it is not held (95) | strong |
| CL8 Values marks | points, log scale, reference lines (130) | Sonnet |
| CL9 Small items | snapping to codons, EC-149 (178); the cursor across block gaps, EC-154 (179, 142); the Antibody view's invalid labels (37); built-in types' fallbacks (143) with the version-aware validate; the views' fixtures through MessagePack (39) | Sonnet |
| CL10 Novice decisions | F3, F9, F11, F14, a test a decision | Sonnet |
| CL11 The drag rule | a drag's work bounded by the screen (the PAE drag's card, the 3,000-contig drag), a perf-tier bound | strong |
| CL12 M5's leftovers | the 200-copy split patch, the Signal zoom frame, Cmd+wheel on the ring at the fit | strong |
| CL13 Carried small items | the BLOSUM62 text; a folder-graph test; a lanes band's history keyed by item; `rowsValue` over a lazy alignment; `seq.render.pairs.filter` | Sonnet |
| CL14 Scaffold checks | `features-planned` over a brief's numbered items; the heavy lane's tight bounds isolated in `vitest.config.ts` | Sonnet |

## Order and parallelism

Three lanes at once (the loop), never two in one file.
- **Wave 1**: CL14, CL1, CL2. CL14 first, so the later briefs are checked by it.
- **Wave 2**: CL3, CL4, CL10. CL10 after CL1: both touch the legend's and the headers' code.
- **GATE 1.**
- **Wave 3**: CL5, CL6, CL9.
- **Wave 4**: CL7, CL8, CL13. CL8 after CL6: both touch the tracks and their layers.
- **Wave 5**: CL11, CL12, apart, so their timings are read on a quiet machine.
- **GATE 2.**
- Then the use gate, the novice pass, the measurements, the last gate and Finalize.

A lane starts when its wave's earlier lane in the same files has merged, not when the whole wave has.

## Stages

- **Stage 0**: the rows above listed and decided; the 14 questions designed, with their mocks (B10); this plan agreed.
- **Then** lanes, gates, the use gate, the novice pass, the measurements and Finalize, as the method's loop says; the
  lanes are written into this plan and [checklist.md](checklist.md) at stage 0.

The items are [checklist.md](checklist.md); this plan holds no status (the method's plan template).

## Re-plans

- 2026-11-23: the 14 design questions take no mocks (the author: "i think we can design/implement these directly. the
  mocks are generally designed for more complex stuff"; then "yes, write it that way and draft the answers"). Each is
  answered in its design and its package's spec; the three that are interactions (a layer's legend, items 5 and 13;
  the stem menu, 8, 26 and 28; the square's way home, 6) are agreed in text, gesture by gesture, before a lane. Items
  8 and 26 are one question. The drafts are [drafts/novice-questions.md](drafts/novice-questions.md).
- 2026-11-23: agreed at GATE 0, with its lanes written in (the author: "Write them in"); the "Draft" line above is
  history. The work folder is `close-out`.
