## Seq* (next): host data

- Seq*: the Storybook's fake structure viewer serves PDB 2OCJ to a Seq* view through a data adapter (`withAdapter`)
  with no file read; its Drop chain edit reopens the source in place, and the view lists the issues a dropped unit
  raises.
- `withAdapter` is generic over a registry's entry kinds, so a viewer's registry passes through; an entry of kind
  `format-runtime` without a runtime is kept as it is.
