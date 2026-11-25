## Seq* (next): an mmCIF runtime

- An opt-in mmCIF runtime (`@molstar/seqx-core/runtimes/mmcif`) reads entities, units, author numbering, links, glycans
  and per-residue values (B-factor, pLDDT) into sequence, molecule and complex nodes, bounded against malformed files;
  its parser loads only on first open.
