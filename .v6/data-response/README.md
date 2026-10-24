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

- An optional param `asResponse`, `true` when unset (`asResponseOf`: `asResponse ?? true`). With it, `apply` downloads (or reads the file) as now, holds
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
- Not a raw `fetch` wrapper in `data-source.ts`: it would bypass the URL cache and ref counts and redo `fetchData` (abort,
  Node path, headers and body). Streaming, headers before the body or byte ranges, if needed, come as a new `DataType`
  through `fetchData` and the manager.

### Consumers

- Every parser takes Response Data, so that `asResponse` can default to `true`. The io parsers are unchanged: the
  transform converts at the top of `apply`, `SO.Data.readText(ctx, a)` for a text format (a `StringLike`: past 512 MB
  a `ChunkedBigString`, as before) and `SO.Data.readBytes(a)` for a binary one. A text parser `from: [String, Binary,
  Response]`: Binary Data is the bytes of a text, decoded, and the `from` contains `String | Binary` (what `RawData`
  and code typed before give), so applying it to such a ref compiles. A binary parser `from: [Binary, Response]`:
  String Data is no binary format; applying one to a ref typed `String | Binary` no longer compiles (a format
  provider's `parse` types its data with Response, so it does).
- `ParseCif` hands a response's bytes to `CIF.parse` (text or BinaryCIF by the first byte), never decoding them;
  `ParseJson` decodes and parses them; `ParsePly` reads the bytes (ASCII or binary by the header); `ParseDx` reads a
  response named `.dxbin` as binary, any other as text (both start with the same text header).
- Binary readers that write into their input or keep views of it read a copy (`readBytes(a, { copy: true })`): DSN6
  flips the byte order of its whole `source.buffer`, TRR swaps it, DCD's frames and VTP's arrays are views of the
  input. CCP4 (through a `FileHandle`), MTZ, XTC, NetCDF, Simularium, Artiatomi EM, G3D, PLY and DX read only.
- `ReadResponse { isBinary }` turns a response into String or Binary Data, for a custom transform that does not read
  one; `DeflateData` takes responses (and copies bytes it passes through).
- One form at a time: `ResponseBody` holds its bytes (shared with the asset manager), or its text once a reader asks
  (`text`, `textNow`, `read(ctx, 'string')`): decoded once, the bytes let go of, by the body and, through `onText`, by
  the asset manager (`letGoOfBytes`, `core/src/util/_spec/assets.md` §6: a file keeps its `File`, a download's bytes
  become a `Blob`). So a text format holds one copy, the text. `bytes` after that encodes the text again and holds
  bytes again (a copy; two readers wanting different forms of one data would flip it, none do). Why: the memory of
  text formats is master's again, CIF's stays one copy of the bytes.
- Sniffing never converts: the body is a `Sniffable` (core `util/sniffable.ts`): its `header` (the first 64 KiB,
  kept in either form) for magic numbers and first lines, `sniffBytes` (the header at an offset), `sniffText` (a
  search over the held form: `asciiIndexOf` on bytes, `includes` on text), `heldBytes` (the bytes when held, for a
  binary format read whole). `DataFormatProvider.applicable(info, o)`, used by `DataFormatRegistry.auto` and Prism's
  `providerFor`, gives a response to `isApplicable` as itself for a provider that sets `sniffsResponse` (the CIF
  ones), as its bytes for one that reads the extension as binary, as its text (turning the body to text: its parser
  reads text anyway) for the others. `guessCifVariant` and `looksLikeMmcifParticles` take a `Sniffable`, a string or
  bytes. Prism sniffs a response's held form (`withContent`, the particle sniff).
- `prism/viewer/src/sources/data-by-url.ts` indexes response cells and returns their held bytes (a copy of a text
  encoded, the body keeping its form).

### Decompression

- `Download` and `ReadFile` take `decompress: 'auto' | 'none'`, unset is `auto`, for text, binary and responses: gzip
  data (told by its bytes) is inflated, cached in the asset manager (`core/src/util/_spec/assets.md` §5). `update`
  recreates on a change. Why: servers host `.gz` files without `Content-Encoding` (the EMDB maps).
- `DeflateData` passes bytes that are not gzip through, so `Download` then `DeflateData` (a pipeline written before)
  still works. The EMDB loader (`packages/extensions/plugin/src/loaders.ts`) no longer applies it; `DownloadFile`'s
  `gzip` format reads the inflated download by its name without `.gz`.
- `DownloadBlob` takes the same `decompress`, and `asResponse` for the whole blob: it gives `Data.ResponseBlob`, each
  entry `{ id, data: ResponseBody }` (one shared hold of each entry's bytes), not a `'response'` kind in `BlobEntry`
  (widening that union would break code that switches over its kinds). `ParseBlob` takes it and parses each entry with
  `CIF.parse` (text or BinaryCIF by the first byte); `DownloadStructure`'s blob path (`asTrajectory`) reads with
  `asResponse` when mmCIF reads responses. `fetchDataMany`'s sources take `shared` and `decompress`.
- `ajaxGet` and `ajaxGetMany` are `fetchData` and `fetchDataMany` (they send a POST with a `body`, so not "get"); the
  old names, and `AjaxGetParams`, `AjaxTask`, `AjaxGetManyEntry`, stay as deprecated aliases.
- Not in scope: zip by its bytes (`PK\x03\x04`), brotli (no native decoder in browsers), bzip2, xz.

## What it costs and saves

- Text CIF by response: one copy of the bytes (held, shared with the parsed file), no decode, no encode, no 512 MB
  limit. Measured on the plugin's path (a file asset, `resolve` and `CIF.parse`, the parsed file and the
  coordinates kept; Node, Apple M4, medians of 3 to 5 runs, retained = heap and external after a full GC):
  - flat-192-cif (124 MB, 1.12M atoms): retained 388 → 270 MB; read 28 → 17 ms, parse 314 → 308 ms, total 353 → 336 ms.
  - 7qpd (1.8 MB): retained 6.3 → 4.6 MB, total 11.5 → 10.2 ms. 1cbs (190 KB): 1.3 → 1.1 MB, 3.1 ms either way.
- The text formats on the byte tokenizer (`TextByteTokenizer`: XYZ, GRO, MOL, SDF, MOL2, PSF, TOP, PRMTOP, LAMMPS,
  cube, DX, OBJ, PLY, Dynamo tbl, CSV) read the bytes (`readBytes`), held once as CIF's. The others (PDB, NDJSON, XSC,
  Kinemage, MVSJ, JSON) hold their text, a download's bytes as a blob: as before.
- A `Download` with `asResponse` that nothing reads still downloads, as today.
- Snapshots saved with `asResponse` load in older v6 builds for formats that accept a string (CIF), not for
  binary-only ones.

## Steps

1. Done. Asset manager: held bytes, `shared` reads, a file decompressed once for them; tests.
2. Done. `Data.Response` / `ResponseBody`, `ReadResponse`, `readText` / `readBytes`; `asResponse` on `Download`,
   `ReadFile` and `DownloadBlob`, the builder overloads.
3. Done. `ParseCif`, `ParseJson`, `ParsePly`; the CIF providers set `sniffsResponse`; `DataFormatProvider.applicable`;
   `data-by-url.ts`; byte sniffing.
4. Done (above).
5. Done: every parser takes a response, `asResponse` defaults to `true` (Consumers); the input-write audit is applied.
6. Open. Move the rules to the specs (the asset manager's are in `core/src/util/_spec/assets.md` already, as the code
   cites it), `CHANGELOG.md` and a `.v6/changelogs/` entry; delete this folder.

Decided: `asResponse` defaults to `true` and every parser takes a response; a custom parser on top of a download
adjusts. The run-time and TypeScript breaks are recorded in `.v6/plans/breaking-v6-changes.md`.

## Where the rules go

- Asset manager (held bytes, shared reads, release): a new `packages/core/src/util/_spec/` (`assets.md`), like
  `core/src/state/_spec/`.
- `Download` / `ReadFile` / `DownloadBlob` `asResponse`, `Data.Response`, `readText` / `readBytes`, `ReadResponse`,
  `sniffsResponse`, `DataFormatProvider.applicable`, the one-form body:
  `packages/plugin/core` has no `_spec`; open: start `packages/plugin/core/_spec/data.md`.
- `ParseCif` reading shared bytes: `.v6/designs/cif-parsing.md` (`CIF.parse` and the plugin), replacing its "the plugin
  keeps reading `.cif` as strings" and the byte-backed `StringLike` pointer; the same pointer in
  `13-lammps-large-text.md`.
