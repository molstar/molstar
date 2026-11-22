# Drift: a merged run of mixed categories drawn neutral, so ClinVar at TP53's fit is mostly grey

**What**: the ranges layer merges a lane's items closer than a pixel into one run (`render/ranges.ts`, about line
508). VC made a run keep its fill when every merged item's fill is equal; a run of mixed fills is drawn in the neutral
range colour. At TP53's fit (19 kb in about 460 px, about 40 bases a pixel) ClinVar's variants merge heavily: 51 of the
band's 85 rects are mixed, and the canvas is about 13,200 px of neutral against 720 of SNV blue
(`seq-next-features-genome-formats-vcf--tp-53-clin-var`).

**Where it disagrees**: genome-formats.md Defaults 2, "a variant a mark at its REF extent, coloured by `category`".

**Options**: (a) a mixed run takes the majority category's fill (ties by legend order), the hover listing the
categories under the pointer; (b) a mixed run drawn as stacked stripes, one per category, by share; (c) variants never
merged below a zoom (more lanes, a taller band); (d) neutral, as now, Defaults 2 amended to "coloured where one
category".

**Recommendation**: (a): the picture reads by category at the fit for the cost of one count a run; (b) is noisier at
1–2 px a run. A fix lane's item in the ranges layer after the author's yes.
