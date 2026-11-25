## Seq* (next): archives

- A view is saved with its data as a zip (`writeArchive`) and opened again (`readArchive`, or `archiveReader` as the
  assets port). Bounded by entries, total and entry size, the compression ratio and names, all read before anything
  is inflated; what is written opens under the same bounds.
