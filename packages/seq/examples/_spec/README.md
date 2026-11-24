# `@molstar/seqx-examples`: the datasets and their views

The example datasets of Seq\* (next), private: the fixtures of each dataset (`fixtures/<dataset>/`) and one builder
module per view (`src/<dataset>/<view>.ts`, exporting `build()`, the view's spec tree written with the language
builder). `views.ts` lists every view as an `ExampleView` (`dataset`, `id`, `title`, `summary`, `build`, and a
`registry` for a view that needs a runtime beyond the default's). The tests of the core and the viewer and the Storybook
import it, so it holds what they share; it is never published.

This spec holds the package's rules, each a rule and its reason; code comments cite its sections.

## Dependencies

- **`@molstar/core`, `@molstar/seqx-core` and `@molstar/seqx-ext-rondo`** (the genome ring's views). Never the viewer: a
  view is a spec, and what draws it is the consumer's. Why: the core's tests import the examples, so a dependency the
  other way would be a cycle.
- **`@molstar/io` is a devDependency**: the record scripts read mmCIF through its text CIF parser (`scripts/cif.ts`),
  and the medium tests that derive from a structure import that script; no view imports it.
- **Entries**: `.` (the views) and a `./*` wildcard onto the modules, so any module is reachable. `.` pulls in nothing
  of Node; `./node-assets` (the Node assets adapter, `node:fs`) and `./trna/generated` (`node:child_process`, `node:fs`)
  are Node-only, so a browser bundle imports neither. `./_test/*` is `null` in the exports: the tests are never
  published.

## Fixtures

- **Small and committed, read through the source protocol like any other data**: a file under 1 MB, a dataset's folder
  under 2 MB (a test checks both). Why: the repository carries them and the Storybook ships them as assets; a dataset
  that needs more is a decision of the author.
- **Real entries, read as published, converted only where no reader exists**: a fixture is in its own format when the
  core reads it (FASTA, Stockholm, Newick, tables, JSON, MessagePack, GFF3, BED, VCF, bigWig, BGZF pages), else in `seq`
  (the protocol's shape), converted at recording (mmCIF, GenBank, a service's metadata). A JSON fixture is trimmed,
  never reshaped, so a spec's paths read the live response too. Why: a reader without Seq\* can take a published file; a
  format the core will not read is no reason to commit megabytes it would never open.
- **Recorded by hand, never in a test**: `scripts/record*.ts` fetch each entry, trim it to what the views read, convert
  it and write the fixture and its provenance. Why: recording is the only step that trusts a remote host, and a test
  never reaches the network. The scripts import the package's rules built (`lib/`), since Node maps no `.js` import to a
  `.ts` file; a run with no dataset only loads the recorders and prints the usage.
- **A real slice of a large indexed file is recorded as pages**: the 64 KiB pages the runtime read, keyed by file and
  page (`<file>.pages.json`), so what the test reads is the container the real tools wrote; a read the recording did not
  make fails naming its page, and a change to a runtime's reads is a re-recording. Why: a trimmed file would test the
  trimming, not the BGZF blocks and the index.
- **Provenance is data**: `fixtures/<dataset>/sources.json` has a row for every file of the folder: the entry and its
  version, the URL, the date recorded, the SHA-256 and size of the download and of the fixture, the licence, the
  attribution, what was trimmed and what is derived (computed in a script or in code, not in the source: contacts, a
  converted numbering). A test checks that every file has a row, that each row's size and checksum are the file's and
  that its licence is allowed. Why: a test can check data; prose cannot.
- **Licences allowed** (`AllowedLicences`, the record script's and the test's one list): `CC0-1.0`, `CC-BY-4.0`, the
  NCBI's and the EBI's public-data terms, UCSC's "freely usable" downloads and `LicenseRef-seqx-generated` (what the
  package makes itself, under the repository's licence). Anything else is refused at recording; one is added by the
  author's decision.
- **Generated, never stored**: every large variant and every computed table (a generated alignment, a genome's links,
  synthetic cytobands) is made in code from a seed, on core's `PCG` (`@molstar/core/data/util/hash-functions`), never a
  private generator. Why: a seed is a few bytes, and one generator makes the same data everywhere.
- **A generated file a view reads by location is committed and compared with its generator**: its row in `sources.json`
  says `LicenseRef-seqx-generated` with no URL, and a test fails when the file is not what the generator writes;
  `SEQ_RECORDS=update` writes it again, reviewed as a diff. A missing file fails: it is written on purpose.
- **Facts are checked against the data**: a dataset has tests that state known facts of its entries (a length, a
  conserved residue, a pair, a feature's span), so a recording that changed the data fails a test, not a picture.

## Reading fixtures

- **A spec names a fixture by its relative location** (`fixtures/p53/tp53-gene.seq`), and the same spec runs in Node and
  in Storybook: Node's assets adapter (`nodeFixtureAssets`) reads it from the package folder, the Storybook's Vite
  plugin serves the folder. Why: a view that differed by host would be tested in one and shown in the other.
- **The Node adapter is confined to `fixtures/`**: no network, no absolute path, nothing above the folder, and no link
  that leads out of it (the file's real path must be under the folder's real path). A refusal is the one the browser's
  policy gives, which the source's opening turns into an issue.
- **Examples load their data from files** (FASTA, Stockholm, a saved service response with its transform, `seq`): as a
  reader of Seq\* would. Inline data is only for a few values or where the inline form is what a view shows.

## Views

- **One module per view, `build()` returning the spec tree**, over a dataset's fixtures; a view that has a second
  instance of its shape (a generated variant, a rotated record) takes the dataset as an optional argument, defaulting to
  the fixture, and runs over the instance that carries its edge case, saying so. Why: an edge case claimed by a view
  that never meets it is not kept.
- **The ring's data** (`genomes/ring`): the assembly by inline lengths, the seeded links (held and as a BEDPE source),
  `syntheticBands(sequences, seed)` for tests and stories without the recorded cytobands. Why: every test and story is
  deterministic with no recording.

## Tests

- **A test per view** (`views.test.ts`, over `ExampleViews`): the spec passes the language's checks; it resolves and
  lays out with no issue but the ones the view is known to report (`Known`, by `dataset/id`: a view whose point is a
  fault says so); every source opens with none either; it renders to SVG settled with no issue. A view that takes long
  to draw (`SlowToRender`) renders in the medium lane and the others are its quick representatives; drawing the same
  spec twice giving the same output is medium for all.
- **A few views keep a draw-list record** (`src/_test/records/*.txt`, committed text compared exactly, written by
  `SEQ_RECORDS=update` and reviewed as a diff); the rest are checked by properties (no issue, bounded draw lists, a hit
  where a mark is), not by records. Why: records are read as diffs, so few and at one size.
- **A story that composes more than one consumer has its twin in `src/_test/`**, composing it as the story does in Node
  on a fake clock. Every asynchronous source a view uses is the synthetic twin of its runtime with an injected clock,
  answering at once in tests (`registry({ delayMs: 0 })`).
- **Tests never reach the network**; the context's assets port is the confined adapter.
