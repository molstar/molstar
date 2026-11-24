# Unreleased

- update documentation for transformation matrix format
- fix the structure cache leaking every entry (molstar#1816): expired entries stayed in memory, each cache hit added a
  duplicate, and the size limit undercounted until the heap ran out; entries are now freed on expiry, evicted until a
  new one fits, and one larger than the cache is not cached
- `cacheMaxSizeInBytes` counts approximate memory (about 2.5x the source size for BinaryCIF, 2.8x for text CIF, plus 2
  MB an entry), not source bytes; the default is 2 GB (was 2 * 1014 MiB); a warning at start when it is over half the
  heap limit
- requests are processed at most `maxConcurrentJobs` at a time (new, default 2; every request that found the queue empty
  ran at once before), and `maxQueueLength` (now a command-line option too) is enforced: more requests get 503 with
  `Retry-After`
- concurrent requests for an entry that is not cached share one read; a model's structure is built, and its custom
  properties attached, once while the entry is cached (was on every request)
- requests whose client is gone (disconnected or timed out) are skipped, or stopped between the entries of a query-many
- remote sources (`sourceMapUrl`): `fetchTimeoutMs` (new, default 30 s) bounds each attempt including the body (was a
  hard 500 ms to the first byte: many RCSB entries, 3J3Q and 7ARV among them, failed); network errors are retried (Node
  reports them as the cause of `fetch failed`); a missing entry answers 404 and a failing source 502
- status codes: invalid parameters and query-many specs answer 400 (were 500 or 404), unmapped sources 404 (was 500),
  other failures 500 (were 404)
- graceful shutdown on SIGTERM/SIGINT and on `shutdownTimeoutMinutes`: stops accepting requests, finishes and sends the
  ones in flight, `health-check` answers 503 meanwhile
- PDBe custom properties: the API responses are cached per entry (a new key per request never hit and grew the cached
  entry)
- `model-server-query` exits when done (the cache timers kept it alive) and reads an input once for all its jobs
- `scripts/profile-memory.mjs` and `scripts/measure-retained.mjs`: memory profile of a running server, and the retained
  size per source byte that calibrates the cache estimate

# 0.9.13

- /surroundingLigands: honor `omit_water=true|false` for REST GET requests (boolean parser previously coerced both to
  `false`)
- /surroundingLigands: stop leaking the asymmetric unit's water chain into the result when `omit_water=true` (water
  residues pulled in via struct_conn covale/metalc edges no longer match every other water in the chain)

# 0.9.12

- add `health-check` endpoint + `healthCheckPath` config prop to report service health

# 0.9.11

# SDF/MOL2 ligand export: fix atom indices when additional atoms are present

# 0.9.10

- /ligand queries: fix atom count reported by SDF/MOL/MOL2 export

# 0.9.9

- /ligand queries: fix behavior for alternate locations
- /ligand queries: handle additional atoms more gracefully
- /ligand queries: better error message for UNL
- /ligand queries: treat deuterium/tritium as hydrogen

# 0.9.8

- fix support for chem_comp_bond and struct_conn categories

# 0.9.7

- add Surrounding Ligands query

# 0.9.6

- optional download parameter

# 0.9.5

- Support molstar_global_model_transform_info category.

# 0.9.4

- bug fix for /ligand queries on metal ions

# 0.9.3

- optional transform parameter

# 0.9.2

- assemblyName in /residueInteraction
- /ligand query
- additional export encoding formats

# 0.9.1

- query-many
- Config overhaul

# 0.9.0

- REST API support.
- Swagger UI support.
- Response schemas.
- Bug fixes.
- Refactored config which can now be provided as a seprate JSON file.

# 0.8.0

- Let's call this an initial version.
