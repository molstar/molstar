## Seq* (next): the change validators are their own entry

- `checkSeqChanges`, `SeqPatchValidator`, `seqBoundsValidator` and `seqPatchValidator` move from
  `@molstar/seqx-core/lang` to `@molstar/seqx-core/lang/changes`, so reading a spec carries no editing: the read
  bundle is 150.6 kB minified (target 160 kB), headless rendering 267.9 kB.
