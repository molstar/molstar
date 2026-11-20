## Seq* (next): data ports

- A host's adapter replaces a format's default (`withAdapter`, one adapter a format); `values.field` reads a source's
  channel (`RegionRequest.field`); a host's edits (`OpenedSource.changes`) reopen the source in place with node identity
  kept, and a selection on a removed node is dropped with one issue.
- `@molstar/seqx-core/io/structure` turns a host-neutral structure description into molecules' nodes, keyed by
  identifiers, bounded and checked (`seq.io.structure.*`).
