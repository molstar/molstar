# Lane CL5: its fix round (the review's findings)

**Model**: strong. **Base**: `seq/cl5` as the lane left it (7 commits on `7f1383e35`); rebase on `packages-seq` first
(FX27 merged since: `viewer/src/input.ts`, `config.ts`, the specs). The brief is [cl5-search.md](cl5-search.md); its
rules, hot files and checks stand. The core `package.json`'s `"./search"` export and `"./search/_test/*": null` are the
lane's local edit, applied by the orchestrator at the merge; keep them local for the tests.

**The orchestrator's checks on `seq/cl5`** (2026-11-24): unit 5,098, medium 677, core perf, `tsc`, lint (5 warnings),
format: exit 0. Size: default 1,563.8 kB, over 1,560 (the lane names 1,565: the orchestrator raises it at the merge, with
its cause, under the author's leave); minimal 1,099.7 of 1,100; headless 978.7 of 980.

**Items**, a commit at each, each with a test failing on its parent:

1. **A regex compile cannot hang**: a counted repeat of an empty body (`A((){999999}){999999}`, 21 characters) copies
   about 10^12 times in `glushkov` (`regex.ts` `positionsOf` counts it 0). Fold a zero-position repeat into `eps`, and
   count the automaton's copies against a bound with its code. Proof: the pattern refused or compiled in under a bound.
2. **No extra minus-strand hits at a chunk boundary (regex)**: `kernel-matcher.ts`'s minus buffer starts at the chunk's
   start, so a leftward-longest match is cut there and reported over the previous chunk's full hit (`A+`, `A+C?`, `GA*`
   with `chunkLetters: 7`). Read the overlap before the start for the minus strand. Proof: the oracle comparison at small
   chunk sizes, every pattern class.
3. **The step budget charges the work done**: a regex step costs `chunks × w`, charged `w` (up to 16× under);
   `.{1,126}.*`, `(.|A)*`, `[ACGT]*A` block 3–6 s in one `find` over 2^20 letters; the exact matcher's Horspool kernel is
   O(n·m) at worst, charged ⌈m/32⌉ a letter (`C`+999×`A` over poly-A, 755 ms a chunk). Charge the real cost and yield
   within the design's "about a second". Proof: the four patterns bounded in steps; a perf-tier case.
4. **Mod+g yields with no hits**: put `search.next` and `search.previous` in `input.ts`'s `YieldWhenIdle`, as
   `view.back`. Proof: a test that Mod+g with no hits is not taken.
5. **`save` of a stale result** refused (`seq.search.save.stale`) or said; the design chooses (Decision 5: "results may
   be out of date"). Proof: a test.
6. **`maxRegexDepth` configurable** (`config.search`, `boundsOf`), in the Configuration table.
7. **The tests**: a cancel while reading and two runs racing; the progress test asserting it is said once after 1 s; the
   EC-159 test asserting the first hit in the view; the batch yield through the viewer's host, not a global `setTimeout`.
8. **The design matches the build**: `design/search.md`'s step commands `session` (with the reason), `Mod+Shift+G`, the
   Configuration rows built, `seq.search.pattern.syntax`, the stale save. Proof: `designs-in-specs` green.

**For the look** (the orchestrator): the p53 and plasmid Search stories, light and dark; `search` and `search-current`
share `token:cat-3`, differing by style only (Defaults' Provisional: settled at the look).

**Followups** (to the TODO at the merge): rows and cells layers not searched; the alignment representative and heavy
bounds not built; `seq.search.save.rects` unused; a source not open yet skipped, not waited for; the Search mock not
built; feature 182 "building".
