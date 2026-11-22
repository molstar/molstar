## Seq* (next): genome formats

- `@molstar/seqx-core/runtimes/genbank`: GenBank flat files read whole, each record a sequence (circular records cyclic) with its features as ranges (`join`, `order`, `complement`, partial ends, a join across a circular record's origin as one item), found by name; Your own file reads `.gb` and `.gbk`.
