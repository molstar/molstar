## Graphics: faster molecular surface

- Molecular surfaces compute their field 1.6–2× faster with the same output: the occlusion tests scan compact
  neighbour lists, pruned per torus, and run only for voxels that would change (7qpd at 0.5 Å: 1.26 s → 0.80 s; 186k
  atoms at 1.3 Å: 1.7–2.0×) (#1988)
