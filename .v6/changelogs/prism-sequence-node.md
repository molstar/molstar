## Prism: each structure's sequence in the structure tree

- Prism: a Sequence node under every structure with polymers, right under Contents (sequence.md §8.1; off with
  `sequence.tree: false`). Its eye picks the structure the band shows (one at a time; Auto, the old follow rule, until the person picks one), its children are the tracks (each an
  eye, its colorings, Remove), its menu has Add track, Chains, Style, Letter colors, Follow automatically and Legend,
  and its inspector the chains (All, None, Polymers, then each entity and its chains). Tracks say what they draw
  ("Letters · Zappo").
- Prism: the chain filter: every chain, polymers only, every copy apart, or the chain instances chosen, written as the
  section's blocks by one patch (the document and its link kept); `sequence.chains` and the `sequence.chains.*`
  commands.
- Prism: tracks and chains are each structure's own (saved state version 4; version 3's hidden tracks are read as the
  followed structure's). `sequence.track.add` / `.remove`, a `structure` on `sequence.track.hide` / `.show`, and
  `sequence.get` lists the chains and the tracks that can be added. Edits made while a structure's document is still
  being made land in it.
- Prism: track providers (`SequenceTrackProvider`, `spec.sequence.tracks`): the way an extension adds sequence data. A
  provider says what its track is and when it applies; the viewer puts it in and takes it out, per structure. The
  ruler, letters and secondary structure are the viewer's own providers; model-archive's pLDDT track is one.
- [Breaking] Prism: `ViewerExtension.sequence` is `{ tracks?, behaviors? }` (was a list of behaviors);
  `SequenceShown.after` is optional.
- Prism: the Sequence button's menu is No sequence, Text, Tracks and Edit in structure tree (`sequence.reveal`): the
  letter colors, Legend and Visible tracks moved to the Sequence node.
- Prism: the band's focus ring shows from the keyboard only, not after a click.
