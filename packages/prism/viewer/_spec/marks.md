# Marks

A view's residue ranges shown in the 3D view and the sequence until cleared (the PAE plot's box; later a contact map, an
alignment). The code is the truth for how (`src/sources/marks.ts`, `api/highlight.ts`, `sequence/link.ts`,
`packages/extensions/model-archive/src/pae/`).

- **View state**: not the selection, the highlight or a style layer; never a history step, never in the state tree, a
  snapshot or a session.
- **Set whole per owner**: `highlight.mark({ owner, marks, dimRest?, focusSequence? })` resolves each target
  ([targets.md](targets.md)), an empty one fails `not-found`; `highlight.unmark({ owner })` clears. Both `ui-only`: the
  agent highlights or selects instead.
- **Two tones, `highlight` and `accent`**: the two styles of a Seq\* highlight channel, so a mark is one color in 3D,
  the sequence and the view (`marks.colors`, per canvas theme). Later marks paint over earlier; no third overlap color
  until Seq\* channels have one.
- **`dimRest`** greys the structure's other residues under the marks, as AlphaFold DB does, so the marks stand out of
  any coloring.
- **`focusSequence`** focuses the sequence on an owner's new marks, all together (Seq\*'s focus zooms to their hull, as
  a selection's: a PAE box's two ranges, often in two chains, both in view). Only new marks move it (by owner, id,
  structure and label: `highlight.mark` makes new objects each call): the same marks set again, or those already there
  when the link is made, are shown where the view is.
- **3D**: an overpaint over each representation's own, put back after Mol\* re-applies overpaint, once per frame,
  repainting only what changed. **Sequence**: a channel per mark (`marks:<owner>:<id>`).
- **Hover**: `highlight.set({ target, sequence: true })` lights residues in the sequence as a 3D hover (`marks.hover`),
  not echoed back to 3D.
