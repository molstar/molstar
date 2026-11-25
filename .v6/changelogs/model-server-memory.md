## ModelServer: memory leak and request handling

- The structure cache freed nothing: expired entries stayed in memory, every hit added a duplicate, and the size limit
  undercounted until the heap ran out (molstar#1816). 35 RCSB entries, 2 cycles, 512 MB cache, after the TTL and a full
  GC: 240 → 416 MB heap and 0.6 → 1 GB outside it retained, growing per cycle → back to the 28 MB baseline. The issue's
  3J3Q and 7ARV, left alone after the TTL: ~300 MB stayed → freed by V8 without a forced GC.
- `cacheMaxSizeInBytes` counts approximate memory (×2.5 BinaryCIF, ×2.8 text CIF, + 2 MB an entry; measured), not
  source bytes.
- Requests run at most `maxConcurrentJobs` (new, default 2) at a time; `maxQueueLength` is enforced (503 +
  `Retry-After`); jobs of clients that left are skipped; concurrent requests for an entry share one read; a model's
  structure and custom properties are built once per cached entry.
- Remote sources: `fetchTimeoutMs` (new, 30 s, body included) replaces a hard 500 ms to the first byte that failed most
  large RCSB entries (3J3Q, 7ARV); network errors are retried; 404 for a missing entry, 502 for a failing source.
- Status codes 400/404/500/502/503 (every failure was 404, bad parameters 500); graceful shutdown on SIGTERM/SIGINT;
  `model-server-query` exits when done.
