# Data.Response: download now, choose the form later

**Temporary.** A working note for branch `v6-data-response`. When the work lands, its rules move to the specs listed in
[Where the rules go](#where-the-rules-go) and this folder is deleted.

Goal: text CIF from a URL or a file is parsed from its bytes, held once, never decoded to a JS string, and nothing that
reads `Data.String` or `Data.Binary` today changes.

## Today

- `Download` and `ReadFile` (`plugin/core/src/state/transforms/data/fetch.ts`) choose the form before reading, by
  `isBinary`: `Data.String` or `Data.Binary`.
- A `.cif` arrives as `Data.String`: bytes, then `TextDecoder` to a string (`utf8ReadLong`), then `parseCifText` encodes
  it back (`utf8Encode`) for the byte tokenizer (`.v6/designs/cif-parsing.md`). flat-192-cif: 124 MB string plus 124 MB
  of bytes retained, about 15 ms of encode, the decode, and the 512 MB string limit.
- The asset manager (`core/src/util/assets.ts`) holds a URL's bytes as a `File` (`new File([data], 'raw-data')`, a
  copy), and every `resolve` reads it again into a new `ArrayBuffer` (`readFromFile`, another copy). A file asset holds
  the user's `File`, and a read of `x.cif.gz` decompresses it on every read (`getCompression` by name).
- `ParseJson` decodes to a string, then wraps it in a `Response` to parse.

## Design

### `Download` / `ReadFile`: `asResponse`

- A new optional param `asResponse` (default `false`). With it, `apply` downloads (or reads the file) as now, holds
  the bytes and returns `Data.Response`: the I/O is eager, as `fetch()` starts its request, and only the choice of form
  and its decode wait for the consumer. Why eager: the download's cost, progress and errors belong to the `Download`
  node (a 404 fails `Download`, not `ParseCif`). `ReadFile` reads eagerly too: a `File` is a handle that can fail
  later (`NotReadableError` once it changes on disk, a missing path in Node, a bad `.gz`), and its decompression runs
  once. A reader that wants part of a large file (trajectory frames by `Blob.slice`) gets its own read mode later, not
  a lazy default. Without the param nothing changes. Why a param: no new transform name, so old snapshots, the classic UI's actions and code that
  matches `t === Download` keep working; a build without the param ignores it and downloads a string, which `ParseCif`
  still accepts.
- The output type already depends on params (`isBinary`): `update` adds `asResponse` to its `Recreate` condition.
  `ReadFile`'s `update` does not recreate on `isBinary` today; it recreates on both.
- `isBinary` is ignored with `asResponse` (hidden in the params UI, said in its description).
- `builders.data.download` / `readFile` get an overload: `asResponse: true` returns a `Data.Response` selector, else the
  current `String | Binary` union. Why: the inferred union would otherwise widen for every caller and break external
  code that assigns `.data` to `string | Uint8Array`.

### `Data.Response`

- `PluginStateObject.Data.Response` holds a `Data.ResponseBody`: the held bytes (decompressed, `bytes`) and
  `read(ctx, type)`, `type` a `DataType`, the result typed by `DataResponse<T>`; `text` and `json` are helpers over it.
- `apply` takes one hold, `resolve(asset, 'binary', { shared: true })` (a view of the manager's bytes: no copy, no
  decode), and `dispose` gives it back, as `Download` does with its wrapper today. `read` serves every form from that
  view and takes no reference of its own: `'binary'` is the view, read-only (a consumer that writes copies it),
  `'string'` decodes it (`utf8ReadLong`), `'json'` decodes and parses it (`JSON.parse`, both native). Why one hold: no
  consumer manages refs, and a consumer done with its read cannot drop the bytes under another one.
- Reads repeat, in any form (`'string'` and then `'json'`), unlike a `fetch` `Response` body (single use); the doc
  comment says so. Each read decodes again: only the bytes are held, not the decoded forms.
- Not a raw `fetch` wrapper in `data-source.ts`: it would bypass the URL cache and ref counts and redo `ajaxGet` (abort,
  Node path, headers and body). Streaming, headers before the body or byte ranges, if needed, come as a new `DataType`
  through `ajaxGet` and the manager.

### Consumers

- `ParseCif` takes a response (`from: [String, Binary, Response]`) and hands its bytes to `CIF.parse`, which tells text
  CIF from BinaryCIF by the first byte. `ParseJson` (`[String, Response]`) reads `'json'`; `ParsePly`
  (`[String, Binary, Response]`) the bytes, telling ASCII from binary by the header.
- The other parsers do not take a response. Why: a `from` of `[String, Response]` or `[Binary, Response]` neither
  contains nor is contained in `String | Binary`, so TypeScript refuses to apply the parser to a `String | Binary` ref
  (what `Download`'s output is typed as): source code that compiles today would not. CIF and PLY already list both. A
  text parser would gain nothing anyway (bytes plus the string while the node lives); a binary one a copy for URLs
  only. `ReadResponse { isBinary }` turns a response into `Data.String` / `Data.Binary` for them.
- Binary readers audited for writes to their input (for a later pass): DSN6 flips the byte order of the whole
  `source.buffer` in place, TRR swaps it in place, DCD's frames and VTP's arrays are views of the input; CCP4 (through
  a `FileHandle`), MTZ, XTC, NetCDF, Simularium, Artiatomi EM, PLY and DX read only. DX has no sniff to tell text from
  binary bytes: it would need one before taking a response.
- A provider declares `readsResponse`: the CIF ones (mmCIF, CIF core, segmentation, density server, structure factors,
  the mmCIF particle assembly, Relion STAR). `isApplicable` then gets text as its UTF-8 bytes: `guessCifVariant` and
  `looksLikeMmcifParticles` search them without decoding (`asciiIndexOf`, `asciiStartsWith`, core `utf8.ts`), tested
  equal to the string search.
- `DataFormatRegistry.responseExtensions`: the extensions every provider of which reads a response (`cif`, `bcif`,
  `mmcif`, `mcif`, `star`). `OpenFiles`, `DownloadFile` and `LoadTrajectory` read a file of one (or for a provider that
  reads one) with `asResponse`; `DownloadStructure` for a provider that does. Prism's `load.entry` goes through
  `DownloadStructure`; its `load.url` and `load.file` (`dataCellOf`, which chooses the format again from the content)
  do not yet.
- `prism/viewer/src/sources/data-by-url.ts` indexes response cells and returns their held bytes.

### Asset manager: bytes held once

The rules are in `packages/core/src/util/_spec/assets.md`. In short: a download's entry holds its bytes, not a `File`
of them; a file's bytes are read (and decompressed) on its first `shared` read and held, while other reads of a file
read it each time, as before (holding them would keep in memory what the file keeps on disk);
`resolve(asset, 'binary', { shared: true })` returns the held bytes, read-only; every other type decodes from them
(`decodeData`, `data-source.ts`). Trade-off: a blob can live outside the JS heap (Chrome may page large ones to disk),
held bytes cannot; for data the parser keeps anyway (text CIF holds its bytes), one copy beats a blob and a copy.

## What it costs and saves

- Text CIF by response: one copy of the bytes (held, shared with the parsed file), no decode, no encode, no 512 MB
  limit. Measured on the plugin's path (a file asset, `resolve` and `CIF.parse`, the parsed file and the
  coordinates kept; Node, Apple M4, medians of 3 to 5 runs, retained = heap and external after a full GC):
  - flat-192-cif (124 MB, 1.12M atoms): retained 388 → 270 MB; read 28 → 17 ms, parse 314 → 308 ms, total 353 → 336 ms.
  - 7qpd (1.8 MB): retained 6.3 → 4.6 MB, total 11.5 → 10.2 ms. 1cbs (190 KB): 1.3 → 1.1 MB, 3.1 ms either way.
- Text formats read as strings are unchanged: they are not read as responses (Consumers).
- A `Download` with `asResponse` that nothing reads still downloads, as today.
- Snapshots saved with `asResponse` load in older v6 builds for formats that accept a string (CIF), not for
  binary-only ones.

## Steps

1. Done. Asset manager: held bytes, `shared` reads, a file decompressed once for them; tests.
2. Done. `Data.Response` / `ResponseBody`, `ReadResponse`; `asResponse` on `Download` and `ReadFile`, the builder
   overloads. No `readAs` helper: the parsers that take a response read it themselves.
3. Done. `ParseCif`, `ParseJson`, `ParsePly`; the CIF providers declare `readsResponse`; `OpenFiles`,
   `DownloadFile`, `LoadTrajectory`, `DownloadStructure`; `data-by-url.ts`; byte sniffing.
4. Done (above).
5. Done as far as it goes: the audit is recorded; the single-type parsers stay as they are (Consumers).
6. Open. Move the rules to the specs (the asset manager's are in `core/src/util/_spec/assets.md` already, as the code
   cites it), `CHANGELOG.md` and a `.v6/changelogs/` entry; delete this folder.

Decided: `Download` and `ReadFile` keep `Data.Response` in `to`. A type of their output written out as `String | Binary`
no longer accepts it, a TypeScript break and not a run-time one (`to` is not checked at run time), recorded in
`.v6/plans/breaking-v6-changes.md`; the builders' overloads avoid it for `builders.data.*`.

## Where the rules go

- Asset manager (held bytes, shared reads, release): a new `packages/core/src/util/_spec/` (`assets.md`), like
  `core/src/state/_spec/`.
- `Download` / `ReadFile` `asResponse`, `Data.Response`, `ReadResponse`, `readsResponse`, `responseExtensions`:
  `packages/plugin/core` has no `_spec`; open: start `packages/plugin/core/_spec/data.md`.
- `ParseCif` reading shared bytes: `.v6/designs/cif-parsing.md` (`CIF.parse` and the plugin), replacing its "the plugin
  keeps reading `.cif` as strings" and the byte-backed `StringLike` pointer; the same pointer in
  `13-lammps-large-text.md`.
