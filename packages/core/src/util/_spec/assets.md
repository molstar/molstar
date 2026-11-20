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

## 7. Sessions

- A session (`.molx`, `PluginStateSnapshotManager.serialize`) holds every held asset's file, at most
  `PluginConfig.State.MaxSessionAssetBytes` (1 GiB; a host sets its own, `serialize`'s `maxAssetBytes` a call's)
  together; a save past it fails, naming the largest, before any is read. Why: the zip is made in memory (each file read
  whole and deflated into a buffer of its size, one buffer for the archive, no zip64), and a session that opened without
  a file it needs would fail later, worse. A state (`.molj`) holds none.
- A file that compresses little (a particle pack, compressed data, images, video: by its name) is stored as it is.

## 8. Network options of a read

- `ResolveOptions` carries what a read of a URL needs to be safe in an app that composes the manager: `credentials`
  (`'omit' | 'same-origin' | 'include'`), `redirect` (`'follow' | 'error'`), `signal` and `maxBytes`; `fetchData`'s
  params take the same. A manager **honors each or refuses the read**: it never reads as if the option were not given.
  Why: an app's policy (a lab server behind a login needs `'include'`, a public viewer wants `'omit'`) cannot be fixed
  in the core, and an option a manager silently drops is a hole the policy believed closed.
- `DefaultAssetManager` passes an option to `fetch` only when it is given: **no option, the same fetch as before** for
  every app (`fetch`'s defaults: `'same-origin'`, `'follow'`, no bound).
- `maxBytes` counts bytes as received (after the transport's content decoding) and is enforced while the body streams: a
  declared `Content-Length` over it is refused before the body is read, an undeclared body is cancelled at the chunk
  that passes it. Both fail with `SizeLimitExceeded`. Held bytes (a cache hit, a file's size) over the bound are refused
  too. Why: a check after the read has already spent the memory it guards.
- `signal` rejects the read as aborted and stops the download. A waiter on a shared download leaves alone on its own
  signal (the download goes on for the others); when a downloader is aborted by its own signal, a waiter starts the
  download again (§2).
- Held data is served to the reads it was read for: an asset held from a read with credentials `'include'` is not served
  to a read of the same asset with `'omit'` (found through `tryFindUrl`, say): that read **is refused**, not fetched
  again, since the holder's reference would otherwise be shared across two answers. A read with the same options (none,
  both) is served as before. Held bytes, a cache hit, and a file are also refused over the bound of the read (by their
  size, the source file's when the bytes were let go of, §6).
- `maxBytes` that is not a size (`NaN`, negative) refuses the read: a bound cannot be honored when it is not one.
- `file://` reads in Node.js: `credentials` and `redirect` mean nothing, the signal is checked before and after the
  read, and the file is `stat`ed against `maxBytes` before it is read.
- The URL cache key includes the credentials, redirect rule and bound of the read: reads with others do not share a
  download. Why: another credential can answer differently, and a bound is a refusal that must not become another read's
  failure. So a bounded and an unbounded read of one URL download separately (accepted: bounded reads are the policy's,
  unbounded ones are the app's). Reads with the same options share as before. The cache behavior is otherwise unchanged;
  range reads (§10) key their pages the same way.
- `SizeLimitExceeded` (`util/data-source`, re-exported from `util/assets`) names its `source`, a URL or a file name.

## 9. The conformance suite

- `checkAssetManager(create)` (`assets-conformance.ts`, `@molstar/core/util/assets-conformance`) returns how a manager
  breaks the interface: a read gives data; references are counted and released (also after a failed read); concurrent
  reads of a URL share one download; a failed download is not cached; gzip is inflated on request; the four options of
  §8 honored or the read refused; and, for a manager with `openRange`, range reads (§10: the bytes asked, by page; a
  changed `ETag` fails; a server that ignores `Range` read once, under `maxBytes`). Data and a function, no test
  framework, as `checkCodec`: an app runs it on its own manager in its tests, `DefaultAssetManager` runs it in core's
  (`_test/assets-conformance.test.ts`, with a derived manager that adds credentials, and one that drops the options,
  which the suite must catch).
- It runs on a fake `fetch` installed on `globalThis` for the run (restored after, one run at a time): nothing reaches
  the network, and a manager is observed by the requests it makes. A manager that does not download through `fetch` is
  outside its reach.

## 10. Range reads

- `openRanges(manager, asset, options)` (`assets-range.ts`, re-exported by `assets.ts`) opens an `AssetRangeReader`: its
  `size`, and `read(ranges, signal)` giving one byte array a range, clipped at the end. It uses the manager's optional
  `openRange`; a manager without one (a custom one) is served by resolving the asset whole (`binary`, `shared`, under
  `maxBytes`) and giving views of its bytes. Why optional: a plugin's custom manager keeps compiling and stays correct,
  only not lazy. Why here: indexed formats (bigWig) read a header, an index and a few blocks of a file of hundreds of
  MB.
- `DefaultAssetManager.openRange` reads by page (64 KiB by default, a power of two from 4 KiB to 4 MiB; the author,
  2026-11-20) and keeps pages in a cache by asset and page, least recently used going first past `maxPages` (256).
  Readers of one key share its pages and its requests in flight; the last reader's `dispose` lets them go, and the
  manager's `clear`/`dispose` closes every reader (its reads reject). The key is the URL cache's (`urlCacheKey`: URL,
  headers, credentials, redirect, bound, §8) with the page size; a file's is its asset. Range readers do not hold a
  manager entry (`assets`, a session's files, do not list them).
- A read gathers the missing pages of all its ranges; adjacent ones are one request, runs apart go in parallel (four at
  a time). A range inside one page is a view of it, not a copy; a page of a run of several is cached as a copy of its
  own bytes, so `maxPages` bounds memory (a view would hold the whole run). A request is shared by the reads that wait
  on it; each leaves on its own signal (rejecting as aborted), and the last to leave stops it. A failed or stopped
  request is not kept: the next read fetches again.
- URLs: `Range: bytes=a-b`, the asset's headers, `credentials` and `redirect` as given (§8). Opening reads page 0, which
  gives the size (`Content-Range`'s total, checked against `maxBytes` before the body) and the validator (a strong
  `ETag`, else `Last-Modified`; a weak `W/` tag is none), compared on every later page with the total. No `If-Range`:
  RFC 9110 forbids it with a weak validator, and it is not CORS-safelisted, so a cross-origin read would preflight every
  page (`Range` with one range is safelisted). A `Content-Range` that is not the range asked (start, end) or has no
  total, a total or validator that changed, a 416 on a valid range, a 200 after the first page, or a body longer or
  shorter than its range fail with `AssetRangeError` (by name across copies, as `SizeLimitExceeded`). Why strict (the
  author, 2026-11-20: "fail on change"): a proxy that mangles ranges, or a file replaced under a reader, would otherwise
  give wrong bytes that decode as data.
- A server that ignores `Range` (a 200 on the first page): its body is read once, under `maxBytes`, and held whole by
  the reader's store; every range is a view of it (never a download a range). It is not put in the URL cache.
- Cross-origin: a server must expose `Content-Range` (and `ETag` or `Last-Modified`, to be compared) through
  `Access-Control-Expose-Headers`. One that hides `Content-Range` (a 206 without it) cannot have its ranges checked: the
  asset is read whole once, under `maxBytes`, as from a server that ignores `Range`. One that hides the validator is
  checked by the total alone.
- A request with a body (a POST) is not ranged: read whole, as `resolve` reads it.
- Files: a file whose bytes the manager holds (a `shared` read) is served from them; otherwise `Blob.slice` a run of
  pages, never the whole file. `maxBytes` is checked against the file's size.
- `file://` URLs in Node.js: `stat` against `maxBytes` before any byte, then reads at an offset with the module
  `setFSModule` set (held in `internal/fs.ts`).
- The implementation is in `util/internal/` (`range-pages.ts`: the page stores and the whole-asset fallback; `fs.ts`:
  the module `setFSModule` sets), which the package does not export (`"./util/internal/*": null`): only the names above
  are public. Why: `DefaultAssetManager` needs the stores, and core's `./*` would otherwise publish them.
