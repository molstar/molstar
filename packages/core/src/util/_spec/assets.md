# Assets: the specification

`@molstar/core/util/assets` (`assets.ts`): the asset manager, which holds downloads and files and their bytes,
ref-counted. `DefaultAssetManager` is the implementation; a plugin takes a custom one from `PluginSpec.assetManager`.
The code is the truth for how; this file records the rules the code relies on, each with its reason.

## 1. References

- `resolve(asset, type, options)` takes a reference, and the returned wrapper's `dispose` (`release`) gives it back. An
  entry lives while it has a reference (a static one, from `set(..., { isStatic })`, outlives its last). Why: data
  cells, snapshots and behaviours share assets; a count lets each let go without knowing the others.
- A failed read gives its reference back (perf-audit-2 01/B12): a file that fails to read is not held.
- `store: false` (or `resolve(asset, type, false)`) bypasses the manager: no reference, no cache.

## 2. The URL byte cache

- URLs read through one cache keyed by method, URL and headers (`urlCacheKey`), not by the response type: a text and a
  binary read of a URL share one download. A request with a body is not cached (it may answer differently).
- Concurrent resolves share the download in flight. A failed download is not cached; its waiters get its error. An
  aborted downloader does not fail its waiters: one of them takes over.

## 3. Bytes held once

- A download's entry holds its bytes as received (a `Uint8Array`), not a `File` of them. Why: a `File` made from the
  bytes is a copy, and every read copied it again into a new `ArrayBuffer`.
- A file's entry holds the file. Its bytes are read, and decompressed (`.gz`, a single-entry `.zip`), on its first
  `shared` read, and held until its last reference. Other reads of a file are read from it each time, as before. Why:
  holding the bytes for every read would keep in memory what the file keeps on disk; a reader that shares them (Response
  Data) keeps them anyway. A static file lets go of its bytes after its last reference.
- `'zip'` reads of a file read the file: its held bytes are one decompressed entry, not the archive.
- `AssetManager.Entry.file` stays a `File`: for a download it is made from the bytes on access (snapshots read it).

## 4. Shared reads

- `resolve(asset, 'binary', { shared: true })` returns the held bytes themselves; the caller must not change them, or
  every later read sees the change. Other `'binary'` reads return a copy. Why opt-in: readers that write into their
  input (DSN6 flips byte order in place, TRR swaps it) or keep views of it in their result (DCD, VTP) need a copy.
- Every other type decodes from the held bytes (`decodeData`, `data-source.ts`): a `'string'` read is one native decode,
  with no copy of the bytes before it. A binary read of bytes read or inflated for it alone (a file without `shared`)
  returns them, not a copy.
- A custom `AssetManager` may ignore `shared` and return a copy.

## 5. Decompression

- `resolve(asset, type, { decompress: true })` inflates gzip bytes before the read; not for `'zip'` reads. Gzip is told
  by its bytes (`1f 8b 08`, `isGzip`), not the name. Why: a server may send `x.cif.gz` without `Content-Encoding` (the
  bytes are gzip) or with it (the browser inflated them already): the name says neither. Text never starts with `0x1f`,
  and the binary formats' headers do not start with these bytes.
- Held bytes inflate once, whoever asks first: a download's inflated bytes are kept on its URL cache entry, beside the
  compressed ones, for every asset of the URL; a file's with its held bytes. A failed inflate is not kept. Reads that
  hold nothing (a file without `shared`, `store: false`) inflate each time.
- A gzip that declares (ISIZE) more than its data can hold (deflate shrinks at most ~1032×) is refused before its buffer
  is allocated. The inflate never writes past the declared size, so a stream over 4 GiB (ISIZE is modulo 2^32) fails, as
  before.
- The inflate of held bytes runs in the first reader's task; when that task is aborted, another reader takes over, as a
  waiter of an aborted download does (§2).
- A custom manager may not inflate (`decompresses` unset): `AssetManager.resolveData`, which `Download`, `ReadFile` and
  `fetchDataMany` use, then asks it for the bytes and inflates them itself, uncached.
- `decompress` off (the default of `resolve`) gives the bytes as received. `Download` and `ReadFile` turn it on unless
  their `decompress` param is `'none'`.

## 6. Letting go of bytes

- `letGoOfBytes(asset)` (optional on `AssetManager`, a hint): the asset's reader holds its data in another form (a
  `ResponseBody` turned to text), so the entry holds no bytes in memory. A file keeps its `File` and drops its held
  bytes; a download's bytes become a `Blob` (`File`), its source from then on, and it leaves the URL cache (the next
  download of the URL fetches again; the entry forgets its key, so its release does not reach a new one). Why: a text
  format then holds one copy, its text, as before held bytes; a session still saves the asset (`file`).
- Not while other assets of the URL share the bytes (the URL cache entry has other holds): they read them.
