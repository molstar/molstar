## Core: range reads of assets

- `openRanges`, `AssetRangeReader`, `AssetRangeError` and an optional `AssetManager.openRange`: byte ranges of a URL or a
  file in 64 KiB pages through a shared, ref-counted cache, failing when the asset changes between pages; a server that
  ignores ranges is read once, under `maxBytes`.
- Seq* (next): an opt-in bigWig runtime (`@molstar/seqx-core/runtimes/bigwig`), values by region and zoom level, bounded
  against malformed files.
