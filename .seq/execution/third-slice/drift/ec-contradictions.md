# Drift: edge cases the code contradicts (EO, 2026-11-19)

Found by EO (`ec-outcomes-draft.md`). For the author at the next gate, each with a recommendation.

1. **EC-92** says TSV and CSV cells are trimmed; the records reader keeps text cells as written (`01` stays `01`, ` x `
   stays), numbers parse past spaces (`core/src/data/_test/records.test.ts` "keeps text cells as written"), as EC-49
   wants. Recommend: amend EC-92 to the code (a cell is data as written; trimming loses meaning).
2. **EC-138** says no hover emphasis for external hovers; the code projects a structure's hover onto every section
   (`core/src/render/_test/emphasis.test.ts`). Recommend: amend EC-138 to the code (a host's hover is what linking is
   for), and test the ≥3:1 hover outline it also asks for.
3. **EC-63** says pairs of several kinds at one cell show glyphs, "+" when they do not fit; the code slices the dots at
   one cell (`core/src/render/pairs.ts` `sideBySide`), untested. Recommend: amend to slices (built, readable), with a
   test of a hit per slice; glyphs to `TODO.md`.
