# Range reads in `@molstar/core` (draft for Gate 0, item 4)

Status: a draft for the author's acceptance before lane BW starts. Sources: BW's brief, version-1.md (open item "range
reads"), base-library.md, ports.md (assets row), protocol.md, performance.md (bigWig region under 50 ms), core's
`util/assets.ts` and `util/_spec/assets.md` §8, Seq*'s `context/assets.ts` (`policyAssets`).

## What BW needs

- A byte range `[offset, offset+length)` of a URL asset or a local `File`/`Blob`; several ranges in one call (bigWig:
  header, then the zoom-level index, then the data blocks of a region, which are contiguous runs).
- Caching by block (page) per asset: a pan re-reads the index and neighbouring blocks, so pages are kept, shared by
  readers, and let go when the last holder disposes (ref-counted like the manager's entries).
- Cancellation by `AbortSignal` (a pan supersedes a read).
- The asset policy's checks unchanged: host, scheme, `credentials`, `redirect`, `maxBytes` (here: of the whole asset's
  size, learned from `Content-Range`/`stat`, and of what one call may fetch).
- A server that ignores `Range` (a 200 with the whole body): detected, the body streamed once, sliced into pages, and
  that outcome remembered for the asset (never one full download per range). `Content-Range` checked: the start, the
  end and the total must match what was asked, else the read fails (a proxy that mangles it must not give wrong bytes).
- `file://` in Node.js: `stat`, then `read` at offset (the module set by `setFSModule`); `maxBytes` honoured as today
  (stat before read).

## Proposed surface

All in `@molstar/core/util/assets` (re-exports from `util/assets.ts`; the implementation in a new `util/assets-range.ts`
so `assets.ts` does not grow). Additive only.

```ts
export interface AssetRange { readonly offset: number; readonly length: number }

export interface AssetRangeReader {
  /** The asset's size in bytes, known after open (`Content-Range` total, `stat`, `Blob.size`). */
  readonly size: number;
  /** One result per range, in order, each a copy-free view of cached pages where one page covers it. Ranges past
   *  the end are clipped; an empty result is a zero-length array. Rejects as aborted with the signal's reason. */
  read(ranges: readonly AssetRange[], signal?: AbortSignal): Promise<Uint8Array[]>;
  /** Gives back the reference: pages are let go with the last holder. */
  dispose(): void;
}

export interface RangeReadOptions
  extends Pick<AssetManager.ResolveOptions, 'credentials' | 'redirect' | 'signal' | 'maxBytes'> {
  /** The page size, bytes (default 64 KiB; a power of two, 4 KiB to 4 MiB). */
  pageSize?: number;
  /** Pages held at most for this asset (default 256: 16 MiB); least recently used go first. */
  maxPages?: number;
}

interface AssetManager {
  /** Optional, like `letGoOfBytes`: a manager without it is served by `openRanges` below. */
  readonly openRange?: (asset: Asset, options?: RangeReadOptions) => Task<AssetRangeReader>;
}

/** Uses `manager.openRange` when there is one; else resolves the whole asset (`binary`, `shared`, bounded by
 *  `maxBytes`) and serves ranges as subarrays: correct, not lazy. */
export function openRanges(manager: AssetManager, asset: Asset, options?: RangeReadOptions): Task<AssetRangeReader>;
```

`DefaultAssetManager.openRange` behaviour:

- `Asset.File` (and a `Blob`: `Asset.File` already wraps `_File`): `file.slice(a, b).arrayBuffer()` per missing page; a
  file already held as bytes by `set` is served from them (no second copy).
- `Asset.Url`: pages are fetched with `Range: bytes=a-b` through `fetchData`'s params (needs `headers`, already
  there; `fetchData` gets an internal `range` result path returning status and `Content-Range`, not public). First read
  fetches page 0 (or the requested pages) which also yields the size; adjacent missing pages merge into one request;
  several non-adjacent runs go as parallel requests (concurrency 4).
- 200 instead of 206: the body is streamed under `maxBytes` (as today), stored as the asset's bytes through the
  existing URL cache, and ranges are served from them; the reader remembers `ranged = false`.
- 206 with a `Content-Range` that is not the asked one, a missing total, 416 on a valid range, or a changed total or
  `ETag`/`Last-Modified` between pages (`If-Range` is sent): fails with `AssetRangeError` (a named error, like
  `SizeLimitExceeded`; matched by name across copies).
- The page cache key is the URL cache's key (`urlCacheKey`: URL, headers, credentials, redirect rule, bound), so
  different options do not share. References are counted with the manager's entries: `dispose()` and `release` free
  pages; `clear()`/`dispose()` of the manager invalidate readers (reads reject).

New exports: `AssetRange`, `AssetRangeReader`, `RangeReadOptions`, `openRanges`, `AssetRangeError`. No existing
signature changes; `AssetManager.openRange` is optional, so custom managers (a plugin's) keep compiling. The
conformance suite (`checkAssetManager`) gains a range section run only when `openRange` exists.

## Stays in Seq*

- The policy (`policyAssets`): checks the host/scheme/`root` on the URL the read resolves to, then calls
  `openRanges(manager, asset, { credentials, redirect, maxBytes, signal })`, where the timeout is the signal's. It
  adds `AssetReader.openRange?(uri, signal)` (optional, as `readText`/`readBytes` are) returning the reader; a bare
  `{ reader }` host (a test or a CLI) gets the whole-asset fallback.
- bigWig's format knowledge: header and index parsing, zoom-level choice, block inflate and decode, mapping a region
  at a level of detail (a bin width) to a zoom level, its own decoded-block cache (small, per source).
- `AssetRefused` mapping and issues on the source (a range error becomes an issue, never a throw).

## Tests (quick lane, no network)

A fake `fetch` (`vi.stubGlobal`) serving a buffer with a toggle for: honest 206; a server that ignores `Range` (200);
a wrong `Content-Range` (start, end, total); a changed total between pages; 416; a body longer than `maxBytes`.
Cases: ranges across pages; many ranges in one call batched into the fewest requests (assert the request log);
repeat reads hit the cache (no request); concurrent reads share a page in flight; abort mid-read rejects as aborted
and the next read refetches (a failed page is not cached); `dispose` frees pages with the last holder, two readers hold
them; `File`/`Blob` through `slice` (assert no whole read); `file://` through a fake `setFSModule` (stat checked against
`maxBytes`, offset reads); a custom manager gets the fallback; credentials/redirect options reach the request; the
conformance suite's new section on `DefaultAssetManager`. Each fails on the parent (the export is missing).

## Alternatives

1. **A free-standing `RangeReader` in `@molstar/core/util`** over `fetch`/`File`, not tied to `AssetManager`: no
   change to the manager's interface, least risk (about 1 day less), but a second download path outside the shared
   cache and refs, the hosts' manager (a Mol* plugin's) not consulted, so the policy checks and the options must be
   duplicated. Breaks the ports.md decision ("the asset manager is reused").
2. **Seq*-only**: bigWig reads whole files through `readBytes` (or the runtime does its own `fetch` with `Range`).
   Zero core change, but whole-file reads of a 100 MB+ bigWig defeat the point, and a private `fetch` bypasses the
   policy and the manager. Fine for the fixtures, not for real tracks.

## Recommendation

The proposed surface (an optional `openRange` on `AssetManager`, `openRanges` as the fallback-aware entry point, the
implementation in its own file). It follows the author's accepted decision, breaks nothing, and keeps format
knowledge in Seq*. Questions for the acceptance: (a) the 64 KiB default page (bigWig blocks are commonly 10-100 kB);
(b) failing on a changed `ETag` between pages (strict, recommended) versus ignoring it.

## Decided by the author

2026-11-20: the surface "Accept"; the details "64 KiB, fail on change" (64 KiB pages by default; a changed `ETag` or
`Last-Modified` between pages fails with `AssetRangeError`). BW builds it as drafted.
