## Seq* (next): the structure link's stories

- Seq* storybook: a fake structure viewer (chain table, MolViewSpec state, own grammar, latency) with the host's adapter
  over `structureFacts`, and the stories Hover, Selection, Focus, Chain maps and Color by this track, each with a
  headless twin.
- `StructureSide.onHover`'s callback takes the event's time (`at`); `linkStructure` drops a side hover older than the
  pointer's last own hover, so a late 3D hover never takes the pointer's.
