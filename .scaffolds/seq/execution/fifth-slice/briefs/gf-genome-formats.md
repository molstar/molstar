# Lane GF: the genome formats (stage 0, design only)

**Task**: draft `drafts/genome-formats.md`: each format a runtime an application adds (sources.md 39), an opt-in entry
like `runtimes/mmcif` and `runtimes/bigwig`: BED (`base: 0`, BED12 blocks as pieces; files.md 82), bigBed, 2bit,
FASTA with its index by region, GFF3 and GTF (the `Parent` hierarchy assembled into transcripts: files.md's tension,
346-347), VCF with tabix by region, PAF, chain files. **BGZF and tabix in `@molstar/io`** (Gate 0, 2): their public
surface agreed here (base-library.md 47: what core and io already hold; range reads from the fourth slice's BW). Each
format's bounds (security.md), its conformance suite (one generated fixture, one recorded real slice, Gate 0, 7: each
fetch asked first), its size row. BAM and CRAM are out of scope (plan, Out of scope).

**Reads**: files.md, sources.md, base-library.md, security.md, protocol.md, performance.md, the bigWig and mmCIF
runtimes as built (`packages/seq/core/src/runtimes`), `@molstar/io`'s readers and core's range reads; precedents:
JBrowse's and igv.js's readers, htslib's specifications (BGZF, tabix, VCF), UCSC's bigBed and 2bit specifications.

**Shape**: [design-drafts.md](../../third-slice/briefs/design-drafts.md) (read first: the principles, the self-check, the interaction matrix, precedents, Defaults with a picture, the Configuration table with whose each choice is, edge cases by id, gaps as yes or no questions with a recommendation). Gate 0 is decided ([plan.md](../plan.md), Gate 0): build on its answers.

**Model**: strong (opus), design only. **Where**: the main checkout, read only except your one draft file; no code, no commits. Other draft lanes write beside you in `drafts/`: touch only your file.

**Report**: under 350 words: the decision in brief, the public surface and language it adds, what it amends in which agreed design, the gaps with recommendations, the self-check.
