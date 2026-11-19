# Lane RW: the fix round (sent 2026-11-19)

Branch `seq/rw-rows`, worktree `.claude/worktrees/agent-acb0f5683e671f189`; rebase onto `packages-seq` at `2f0e4e417`
(SY's merge `6b36ad68d` plus the `./story` export). Each fix lands with a test shown failing on its parent.

**The bands seam with SY (decided, within stories.md and grids-and-rows.md)**: the band record's `bandTop` is in
display rows (fractional), the layout turns it into px; history and kept views store the top row by key plus
`rowScale`, restored clamped (a folded key to its summary row, a hidden one to the nearest shown); a fold above the
window keeps its top row. `_viewports` takes the ViewModel's bands (RW, with folds), a shown story frame's
`bandTop`/`rowScale` merged over them keeping folds (SY's TODO names the seam: close it); test a step's motion over a
folded band keeps the folds. Bands get their own snapshot field rather than the `'\u0001bands'` pseudo-group.

**Review findings**:
2. a pinned row's tree leaf, hover band and selection tint use its band index: an in-place display index;
3. image tiles with a pinned band and scroll squeeze rows: pinned rows drawn on their own, tiles offset from `pinned`;
4. `paintRows`'s huge-selection skip unreachable: clamp to rows in view plus pinned (test at 200k lazy rows);
6. aligned FASTA takes the width from the first record: the most common width, the rest reported; duplicate labels
   reported; tests for the non-letter check, `~`→`-`, a row named twice in `matchLeaves`;
7. config to the design's table: `rows.reduce` built or dropped; `rows.rowHeight`, `rows.treeHitPx` fixed; Min/Max
   RowScale named in the _spec; logo `rows: shown` and `summary: 'logo'` built or removed; `numbers: gutters` draws the
   last number; `reveal.opensClades`, `rows.clades.cutLine` built if small, else TODO;
8. `layoutTree` and `measureEntries` allocate per pan: kept on the lay; a tree-order allocation case;
9. hits allocate per hover (scratch, measured in `rows.quick.ts`); the tree hit is the nearest branch within 6 px;
10. `prepareColors` merges the rows' measures over `theme.measures`;
11. the text figure on `layout === 'text'` only; `rowScale` from `display.count`;
12. `input.ts`: a wheel over a `row` hit from a non-rows layer still pans;
14. tests where mutations survived (pinInPlace false, cladesWithin's parent check, a pinned row in a fold, identity's
    gap exclusion, rowAtY in the pinned band, field order, summary alpha, cell-run merge, `summary: 'bar'`, pin toggle,
    hide's merge, the zoom's anchor; viewer tests for pin/hide/show/order/collapseTo; matrix cells as binding tests);
15. one `RangeLabel`/`firstNumberOf`, one `Alphabet`/`ColorFallback`, leftovers dropped;
5. the pyramid: measured tolerance kept (gaps.md 5, the author's); the share checked where the letter is wrong too.

**The look** (light and dark): A. Scale: the SVG export has no `<image>`, so rows are blank at the fit (export the
tiles; test renderSvg at the fit of 10k rows); B. labels cut to 3 characters: the header is the tree's 120 px plus the
labels' column, a label that doesn't fit elided with an ellipsis; C. Folded: wedges cover summary labels (after B).

**Design items not built**: build Shift+click extend, the summary row's double-click/Enter, `seq.render.rows.second`,
pins and hidden on lazy alignments (or an issue), the summary's partial mark; a lazy band below a pixel gets coarse
levels or a neutral fill; the rest listed in `packages/seq/TODO.md` with their sections.

**Changelog text** (from RW's first report, to write as `.v6/changelogs/seq-rows.md` at merge): multiple sequence
alignments: aligned FASTA, Stockholm and Newick readers; `alignment` and `tree` kinds; a `rows` layer (tree, clades
folded into consensus summary rows, pinned and hidden rows, order by tree, data or a field, row zoom `]`/`[`, wheel
scroll, held or lazy rows to 10,000,000, image tiles from a shared reduction); a sequence `logo` layer; rows in the
text figure; Globins examples and stories. No hot-file changes needed.
