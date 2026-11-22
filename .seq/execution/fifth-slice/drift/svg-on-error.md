# Drift: the SVG export blank on a spec with an error, the canvas drawing what it can

**What**: a spec whose `check` finds an error (Features › Builder › Check: a layer's `data` naming `/data/helo`, which
does not exist) draws on the canvas (the letters, an empty hits track) while `exportSvg` returns a blank picture
(`packages/seq/core/src/render/export.ts:435`, `:441`: `if (hasErrors(issues)) return blank(true)`).

**Where it disagrees**: the author's rule that the canvas and the SVG of one spec lay out alike (fixes.md, "The canvas
and the SVG of one spec lay out differently"); the export design's "the SVG is the canvas's picture".

**Options**: (a) the export draws what resolves, as the canvas does, its issues listed in the picture's metadata;
(b) keep the refusal (an export is a document, so a broken spec makes none) and have the canvas say so too;
(c) keep both as they are, the story saying it (BS did this).

**Recommendation**: (a): a host exporting a figure with one bad layer gets the rest, as on screen; the issue stays in
the status line and the picture's `<desc>`. A fix lane's item after the author's yes.
