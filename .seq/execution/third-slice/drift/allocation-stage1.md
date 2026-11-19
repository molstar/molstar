# Drift: a layout allocation regression after stage 1

The allocation tier's "p53 on its gene, whole protein" pan measures 70.5 KB at `3241b2eef` (after TH), against its
64 KB bound and the ~40 KB recorded when the bound was set; AX and SL each reported the tier green, TH did not run it.
TX's interface adds ~5 KB more (75.1 KB) and names per-section closures and contexts as the rest. Recommendation:
bisect AX..TH on that scenario, fix the cause (likely TH's channel pass or the theme objects made per layout), and
keep the bound; at the latest at stage 7's measurement.
