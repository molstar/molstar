## Prism: each structure's sequence in the structure tree

- Prism: a Sequence node under every structure with polymers, right under Contents (sequence.md §8.1; off with
  `sequence.tree: false`). Its eye picks the structure the band shows (one at a time; Auto, the old follow rule, until
  the person picks one), its children are the tracks (each an eye, its colorings, Remove), its menu has Add track,
  Chains, Style, Letter colors, Follow automatically and Legend, and its inspector, drawn as a representation's quick
  settings, has Add track (the command menu at its group), Style, Letters and the chains (All, None, Polymers, then each
  entity and its chains; None takes the band away). Tracks say what they draw ("Letters · Zappo").
- Prism: the chain filter: every chain, polymers only, every copy apart, or the chain instances chosen, written as the
  section's blocks by one patch (the document and its link kept); `sequence.chains` and the `sequence.chains.*`
  commands.
- Prism: tracks and chains are each structure's own (saved state version 4; version 3's hidden tracks are read as the
  followed structure's). `sequence.track.add` / `.remove`, a `structure` on `sequence.track.hide` / `.show`, and
  `sequence.get` lists the chains and the tracks that can be added. Edits made while a structure's document is still
  being made land in it.
- Prism: track providers (`SequenceTrackProvider`, `spec.sequence.tracks`): the way an extension adds sequence data. A
  provider says what its track is and when it applies; the viewer puts it in and takes it out, per structure. The ruler,
  letters and secondary structure are the viewer's own providers; model-archive's pLDDT track is one.
- [Breaking] Prism: `ViewerExtension.sequence` is `{ tracks?, behaviors? }` (was a list of behaviors);
  `SequenceShown.after` is optional.
- Prism: the Sequence button's menu is No sequence, Text, Tracks and Edit in structure tree (`sequence.reveal`): the
  letter colors, Legend and Visible tracks moved to the Sequence node.
- Prism: the band's focus ring shows from the keyboard only, not after a click.
- Prism: the structure tree's eye and ⋯ buttons show the app's tooltips, not the system's.
- Prism: Reset view below the sequence band frames the structure as close as Mol*'s own reset allows (its perspective
  margin counted in), not farther back by the band's height on top of it.
- Prism: hiding the sequence band fades it in place (it used to jump to a panel's place, then fade).
- model-archive: `model-archive.fetch-pae` also gets a Model Archive model's PAE from the associated file the model
  names (`ma_entry_associated_files`, a CIF or an entry of a zip), on demand: a cell under the model, as AlphaFold DB's
  PAE. The structure row's entry and the PAE card's button say "Fetch PAE".
- model-archive: a pairwise metric that names a chain the model does not have is read without an error.
- model-archive: the PAE row in the Structure panel (was "Predicted aligned error") is under every structure whose PAE
  can be fetched (Model Archive, AlphaFold DB), not only those that have it ("Not loaded · Model Archive"): its click
  fetches the PAE, then opens the card.
- Prism: a cut-off tree row's tooltip card is the row's height, its text level with the row's and in the row's sizes
  and weights.
- Prism: the sequence's secondary-structure track is the one the 3D view draws (Mol*'s computed assignment: the file's
  where it has one, else DSSP), so a predicted model (Model Archive, AlphaFold DB) without `struct_conf` has one.
