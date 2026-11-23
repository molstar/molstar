# Lane GB: GenBank read whole

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/gb`. Model:
strong (a file is untrusted input). Runtimes and their stories only, never the language. Beside you: RG (core's store,
`data/`, `render/pairs*`, `grid*`, `layer-data.ts`, the matrix layer) and FX5's fix round (navigators, menus).
Shared files (`runtimes/internal/`, `data/protocol.ts`) only grow; touch nothing of RG's.

## The design section (to add to genome-formats.md as "GenBank (`runtimes/genbank`)", the orchestrator's; follow it)

- **What it answers**: a file of one or more records (`LOCUS` … `//`), each **a sequence** named by LOCUS's name (an
  `ACCESSION.VERSION` kept as an alias by the names rule), its length LOCUS's, `circular` a cyclic system, the
  molecule type (`DNA`, `RNA`, `mRNA`, `AA`) its alphabet, ORIGIN's letters (digits and spaces dropped, case kept as
  read); **its features as ranges** on that sequence: the key the `category`, the qualifiers fields (`/label`,
  `/gene`, `/product`, `/note`, `/locus_tag`, any other by its name as text; `/translation` dropped, counted, said
  once), strand from `complement`; `id` `<key>@<start>-<end>` with `#2` … on a repeat in file order (the BED rule);
  the label `/label`, else `/gene`, else `/locus_tag`, else the key.
- **Locations**: `a..b`, a single base `a`, `complement(…)`, `join(…)` and `order(…)` as pieces of one item (GFF3's
  parts: `itemStart`/`itemEnd` the hull), nesting as INSDC writes it (`complement(join(…))`, `join(complement(…),…)`);
  `<` and `>` partial ends kept as fields `partial_start`, `partial_end`; `a.b` (one base unknown within) read as its
  hull, said once. **A `join` over a circular record's origin** (`join(2600..2686,1..120)`) is one item across the
  origin (the cyclic system's wrap, as built for ranges). **Refused per feature, with an issue, never the file**: remote
  references (`J00194.1:100..202`), `a^b` sites, a location past the record's length or a malformed one
  (`seq.source.genbank.location`).
- **Header fields**: DEFINITION as the sequence's description, SOURCE/ORGANISM kept as the record's `organism`; the
  rest (REFERENCE, COMMENT, DBLINK, KEYWORDS) skipped and counted once.
- **Names**: `findNames` over `/label`, `/gene`, `/locus_tag` (and the record's LOCUS name and accession), the map built
  once at open, as GG's.
- **Bounds** (security.md): read whole within `limits.maxBytes` and `maxRows` (features), a record's letters within the
  2^31 − 1 length, records a file `maxRecords` (1,000), qualifiers a feature 256, a qualifier's text 64 KiB, a
  location's pieces 10,000 and its nesting depth 16; over a bound: the file `seq.source.too_large`, a feature dropped
  and said. ORIGIN letters that disagree with LOCUS's length: the letters win, said.
- **Refused, the file**: no `LOCUS` line, a record without `//` at the end of the file (truncated: the records before
  it kept, said), an `AA` record's letters outside the protein alphabet beyond the built rule.
- **Not answered**: EMBL flat files and GenPept-only fields; CONTIG/`WGS` assemblies of other records (refused, said).
- **Configuration**: `maxRecords`, `maxQualifiers`, `maxQualifierBytes`, `maxLocationPieces` runtime options (defaults
  above), copied into core `_spec`'s table; `format: 'genbank'`; the runtime an opt-in entry.

Write this section, as built, into your report as the text the orchestrator puts in genome-formats.md (a hot file).
If building it shows a rule above is wrong, say so with the reason; do not silently differ.

## Task

- `runtimes/genbank` (`GenbankRuntime`, `genbankRuntime`), on `runtimes/internal/` (lines, issues, interval index,
  answers) as BED, GFF and VCF use it; read whole; its own size row under 30 kB in `packages/seq/_measure/sizes.mjs`;
  the viewers' bounds unchanged (reading 240, headless 870, default 1,400, minimal 1,000 kB).
- **Fixtures**: the plasmid's derived features written as GenBank by `examples/scripts/record.ts` (as its GFF3 is: an
  NCBI-shaped record, bla's `join` across the origin, lacZα `complement`), committed. **The real L09137.2 flat file is
  a download the author has not yet approved: do not fetch it.** Build everything on the generated file; leave a
  clearly named, skipped test for the real file and report what you would fetch (URL, size, licence to check).
- **Stories**: "Try pUC19.gb" in Examples › Services › Your own file loads the generated GenBank (`dropped-file.ts`:
  `.gb`/`.gbk` leave the refused list; a story test that every Try button loads a view); a plasmid story reading the
  GenBank file beside the GFF3 one, **the two views equal** (a twin test comparing their draw lists); a misread file;
  `Features/Genome formats/GenBank`; light and dark, canvas and SVG, the headless twin on a fake clock.
- **Tests**: the location grammar (each form above, nesting, partials, the origin `join`, each refusal per feature);
  multi-record files; `findNames`; fuzzed and truncated files without a throw; each bound refused before it allocates;
  quadratic scans checked by counting (VC's review found three: bound every `indexOf` to its line or field);
  10,000 features parsed under 100 ms (heavy, `itWithin`; 1,000 quick). EC ids the design names for read-whole
  formats (EC-40, EC-51, EC-86, EC-92, EC-195, EC-196 at least) in test titles.
- core `_spec` Data section for GenBank and its Configuration rows (the config-keys check).

**Where**: a worktree on `seq/gb`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at
each small green step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`): list them with the exact edit. Anything in
`@molstar/core` or `@molstar/io`, the language, or the network stops the lane and comes back as a question. A file
past 800 lines is a question; past 1500 a failure.

**Checks**: the quick lane (`pnpm test:unit packages/seq/`, under 2 s a touched folder), `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, `pnpm -s -C packages/seq/core size`, `pnpm test:medium packages/seq/`, `pnpm test:heavy
packages/seq/` only when `uptime` load is under 6 (else say so), the `.scaffolds/seq/_checks` (features-planned with `--plan
.scaffolds/seq/execution/fifth-slice/plan.md`). A fix's test fails on its parent.

**Report**: under 350 words: commits; checks with numbers; story ids; the public surface; the design section as built
(for genome-formats.md); the real file's fetch (URL, size, licence); sizes; hot-file needs; the `.v6/changelogs/` text.
