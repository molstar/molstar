# Drift: edge cases the code contradicts (EO, 2026-11-19)

Found by EO (`ec-outcomes-draft.md`). Finalize (CL, 2026-11-20) resolved items 1 (EC-92 changed, as files.md already
decides), 2 (EC-138's projected hover changed, as viewer-models.md's projection decides) and 4 (a folded clade's summary
row unnamed: against grids-and-rows.md's accessibility, to `packages/seq/TODO.md`), recorded in edge-cases.md's
"Outcomes of the third slice". Two remain for the author, each with a recommendation.

3. **EC-63** says pairs of several kinds at one cell show glyphs, "+" when they do not fit, and grids-and-rows.md
   ("Kinds at one cell") says so too; the code slices the dots at one cell (`core/src/render/pairs.ts` `sideBySide`),
   each slice its own hit (ET's `ec-pairs-kinds.test.ts`). Recommend: amend the design to slices (built, readable,
   tested); glyphs to `TODO.md`.
5. **EC-138's hover outline** holds 3:1 against the page, tracks and residue fills in both themes (ET's
   `core/src/render/_test/ec-hover-contrast.test.ts`), not against the saturated data fills (dark `seq-value` 2.18:1;
   `seq-cat-*`, `seq-ramp-*`, `seq-ss-*` about 1.4–3.4:1; light `seq-ramp-4` 1.76). Recommend: the EC asks 3:1 against
   the track (what the outline sits on), not every fill; a fill it fails gets the double outline (light and dark)
   ux-foundations already uses for the selection.
