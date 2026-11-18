## Seq* (next): security bounds

- Specs and data are untrusted input, and what they size is refused before the work starts: `checkSeqSpec`, the patch
  check, whole roots (`viewer.load` rejects with `.issues`), previews and story steps refuse a spec past 64 MB, 256
  levels of nesting (`custom` and `args` included), 64 nested query forms or a column of 2^31 or more values, with
  coded issues `seq.lang.bound.*`; issues kept are capped (10 places per code, 100 per node); answers and `seq` files
  are refused past a column of 2^31−1 values, a block of 2^24 values or too many dynamic nodes; nested sources default
  to 4 levels; `headlessAssetPolicy({ root, hosts })` reads only files under a root and listed https hosts; the SVG
  backend replaces lone surrogates; `renderSvg` takes `bounds`. The bounds are configuration items
  (`SeqCoreConfig.bounds`, `bounds.*` in the viewer).
