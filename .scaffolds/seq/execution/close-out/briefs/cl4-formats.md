# Lane CL4: a canonical transcript, the VCF's sample names, the matrix's names

**Model**: Sonnet (a closed item list with named outputs: three answers of a sentence each, in readers and a card).
**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl4`. Beside
you: CL1 (the legend, `viewer/src/legends.ts`, `mount-legend.ts`, the themes design) and CL2 (the ring,
`packages/seq/ext-rondo/`, the ring's stories). Never edit their files.

**Why**: the close-out plan, wave 2. The fifth slice's novice passes found 40+ TP53 transcripts with none marked
canonical (10), the VCF's 40 sample rows unnamed but on hover (11), and the methylation matrix's samples called "State"
and its chromosome "chr" (look 34). The answers are agreed (the author, 2026-11-23: "Agree as drafted"), and a collapsed
gene keeps every transcript in version 1 (canonical-only is a TODO line, not this lane).

**Read first**:
- `.scaffolds/seq/execution/close-out/drafts/novice-questions.md`, **F** (10, 11);
- `.scaffolds/seq/execution/fifth-slice/novice-findings.md`, items 10, 11 and 34;
- `.scaffolds/seq/design/genome-formats.md`: GFF3 and GTF; VCF; What every runtime shares; Defaults; Configuration;
  `.scaffolds/seq/design/region-reads.md` (the matrix);
- `packages/seq/core/_spec/README.md` (the runtimes, the matrix layer) and `packages/seq/viewer/_spec/README.md`
  (Hover card);
- the code: `core/src/runtimes/{gff-parse,gff-assemble,gff-read,vcf-header,vcf-genotypes,vcf-read}.ts`,
  `core/src/render/{matrix,rows-columns}.ts`, the card providers (`viewer/src/card-providers.ts`, `describe*.ts`).

**Items**, a commit at each:

0. **The checks first**: how the alignments' rows name their rows today (a names column or rail), and whether the
   matrix layer can take the same part; which GFF3 and GTF attributes the fixtures carry (`tag=Ensembl_canonical`,
   `MANE_Select`). A contradiction with the agreed answer stops the lane: report it as a drift note.
1. **A canonical transcript (10)**: the GFF3 and GTF readers keep the `Ensembl_canonical` and `MANE_Select` tags; a
   tagged transcript is drawn first in its gene's stack with its label bold, and its card says "Canonical" (or "MANE
   Select"). With no tag nothing is marked. Proof: reader tests on a fixture with and without the tag (EC-227's
   transcripts); a layout test on the stack's order; a card test.
2. **The VCF's sample names (11)**: the sample rows get a names part from the VCF header's sample columns, the same
   part the alignments' rows use. The card reads "Sample S11 · 0/1 (heterozygous)", not "State S11 / Value 0/1".
   Proof: tests on "Features › Genome formats › VCF › Genotypes" and "Examples › Genomes › Views › Samples" (names drawn;
   the card's text), with their twins.
3. **The matrix's names (look 34)**: the methylation matrix's samples named by the same part, never "State"; its
   chromosome named as the data names it, never "chr"; a missing cell in dark distinct from the ramp's low end. Proof:
   tests on "Examples › Genomes › Views › Methylation", light and dark.
4. **The designs and specs**: the answers into `genome-formats.md` and `region-reads.md`, their Storybook Design pages,
   and the core and viewer specs, in the commits that build them; a TODO line for "a collapsed gene draws only its
   canonical transcript" named in your report (you don't edit the TODO). Proof: `designs-in-specs` green; the pages.

**How**:
- **host neutral**: names are drawn in the view's own part; nothing of the host's chrome;
- **two consumers**: the names part is tested with the alignment's rows and the VCF's samples live in one test;
- **minimal code**: reuse the rows' names part; no second one; module size (§module-size);
- **its performance representative**: the genotype matrix of 40 samples by region keeps its existing bound; names are
  laid out when the view settles, not per frame;
- **lazy parts awaited**: region reads are lazy: await them before a gesture;
- **the agreed answer** is the mock; stories with headless twins on a fake clock.

**Where**: a worktree on `seq/cl4`. Touches `packages/seq/core/src/{runtimes,render}`, the card providers in
`packages/seq/viewer/src`, the genome-formats and region-reads stories and Design pages, the two designs and the specs.
Out of scope: the legend (CL1), the ring and its stories (CL2), the square (CL3). Headers per the root `AGENTS.md`. Hot
files never staged: `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`, the work folder. Name them in the
report (the changelog line too). The language does not change: a change to it stops the lane. Write only in your
worktree and your own temporary space. No network: fixtures are recorded already; a new one is made in code.

**Checks**, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/ cli/seqx-render/`;
- `pnpm test:medium packages/seq/`;
- `pnpm -s -C packages/seq/core size` (read the rows; the runtimes' rows under 30 kB each);
- `npx tsc -b`; `pnpm lint`; `pnpm format:check`;
- the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first: `git add .scaffolds packages/seq`).

Your quick-lane time per touched folder stays under 2 s. Already red at the base: name it, don't fix it.

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report the action and why; don't work around it).

**Report** (under 400 words): commits; each item with its proof; each check with its exit code and test count; changed
expectations and why; what to look at (stories, light and dark); drift, gaps and followups; the changelog text.
