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
  with no copy of the bytes before it.
- A custom `AssetManager` may ignore `shared` and return a copy.
