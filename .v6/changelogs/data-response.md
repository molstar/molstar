## Plugin: downloads and files as Response Data

- [Breaking] `Download`, `ReadFile` and `DownloadBlob` give `Data.Response` (`DownloadBlob`: `Data.ResponseBlob`)
  unless `asResponse: false`: the bytes, held once, read by the parser as text or binary (`data.bytes`,
  `data.text(ctx)`, `data.json(ctx)`, `data.read(ctx, type)`). Every built-in parser takes String, Binary and Response
  Data; a custom parser on top of a download adds `Data.Response` to its `from` and reads through
  `SO.Data.readText` / `SO.Data.readBytes`, or the download sets `asResponse: false`, or `ReadResponse` goes between.
  `builders.data.download` / `readFile` / `downloadBlob` return Response selectors unless `asResponse: false`.
- Text CIF parses from the bytes, never decoded to a string and encoded back: flat-192-cif (124 MB) retained 388 → 270
  MB, read and parse 353 → 336 ms; 7qpd 6.3 → 4.6 MB. No 512 MB string limit for CIF.
- [Breaking] `Download`, `ReadFile` and `DownloadBlob` inflate gzip, told by its bytes and not the name (`decompress`,
  `'auto'` unset, `'none'` for the bytes as received): a server that sends `.gz` without `Content-Encoding` (the EMDB
  maps) gives the data, not gzip bytes. `DeflateData` passes bytes that are not gzip through.
- Asset manager: a download is held as its bytes (not a `File` copy); `resolve(asset, type, { shared, decompress })`;
  `AssetManager.resolveData` inflates for a custom manager that does not.
- `ajaxGet` / `ajaxGetMany` are `fetchData` / `fetchDataMany` (and `FetchDataParams`, `FetchDataTask`,
  `FetchDataManyEntry`); the old names stay as deprecated aliases.
- `ReadFile` recreates when its file changes (it kept the old file's data).
- A response body holds one form, its bytes or its text once a parser reads it: a text format holds one copy. Its
  header and sniffs (`sniffBytes`, `sniffText`) work on either; core `Sniffable`; a provider sets `sniffsResponse` to
  sniff a response as itself (`DataFormatProvider.applicable`).
- [Breaking] Binary parsers take Binary and Response Data, not String Data: applying one to a ref typed
  `String | Binary` no longer compiles.
- Core `asciiIndexOf`, `asciiStartsWith`; zip `isGzip`, `gzipDeclaredSize`.
