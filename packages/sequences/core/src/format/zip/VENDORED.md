# Mol*'s zip, copied

Source: vendored Mol* 5.11.0 (`1ebe24097`), `src/mol-util/zip/` (`bin`, `checksum`, `deflate`, `huffman`, `inflate`,
`util`, `zip`; UZIP.js by Photopea, ported by Mol*, MIT). The files are unchanged except their import paths (`./x` is
`./x.ts`; the two imports from outside the folder point at `shims.ts`) and the local changes below. Biome and
`seqstar:ported-lines` skip the folder, so it stays as Mol*'s.

Local changes, each marked `STORIES-UNTOLD(Zn)` in the code (and worth sending upstream):

- **Z1** (2026-10-14, the review gate's untrusted inputs GI-01, GI-02, GI-08, GI-09), `zip.ts` and `inflate.ts`: the end
  record is looked for only in the last 64 KB (a truncated archive searched below 0 for ever); `unzip` reads the central
  directory whole before inflating anything, takes its sizes, name and CRC (an archive written with a data descriptor,
  as Finder's Compress writes it, has zeros in the local header), and checks optional `UnzipLimits` (an entry count and
  a total declared size) and that no local entry is listed twice; every entry's bytes must be its declared size and
  CRC-32; `_inflate` stops at the end of the buffer it is given instead of growing it, and returns what it inflated. The
  fallback inflater (no `DecompressionStream`) still grows its buffer before the size check.

`shims.ts` stands in for what the files import from outside the folder:

- `mol-task`: `RuntimeContext` (a context that never asks to yield) and `Task` (`Task.create(name, f).run()` calls `f`
  at once);
- `type-helpers`: `NumberArray` and `assertUnreachable`.

Each shimmed type is also a value, so the files' `import { RuntimeContext }` (Mol*'s style) loads under Node's type
stripping.

Only `readSvsx` and `writeSvsx` use it (seqspec:sources#archives); nothing else imports the folder.

When this code moves into Mol*, delete this folder and its shims and import `mol-util/zip/zip`.
