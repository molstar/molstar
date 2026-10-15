# Marks: a view's residue ranges on the scene

A view other than the 3D view and the sequence sometimes needs to point at residues in both for a while: the PAE plot's
box (its aligned and scored ranges), later a contact map's, an alignment's. Marks are that, once, for every such view.
They are not the selection (the person's, in the history's journal), the highlight (the pointer's, gone on the next
move) or a style layer (in the state tree). Code: `packages/prism/viewer/src/sources/marks.ts`, `api/highlight.ts`,
`sequence/link.ts`.

## 1. What a mark is

```ts
interface Mark {
  owner: string; // who set it ('model-archive.pae'): an owner's marks are set and cleared together
  id: string; // 'aligned', 'scored'
  label: string; // "Aligned A 92–190", for a legend or a tooltip
  tone: 'highlight' | 'accent';
  structure: string; // the structure's state ref
  loci: StructureElement.Loci;
  dimRest?: boolean; // grey the rest of the structure while shown
}
```

- **Two tones**: `highlight` (the sequence's yellow) and `accent` (Prism's accent), the two styles a Seq\* highlight
  channel has (`HighlightStyle`), so a mark is the same color in the 3D view, the sequence and the view that set it.
  `MarkTones` gives them per canvas theme; a view reads them from `marks.colors` (which follows the background).
- **Later over earlier**: where marks overlap, the one set later shows, in the 3D view and in the sequence. A view
  orders its marks so the overlap reads (the PAE plot sets the scored range last: it is painted over the aligned, and
  the plot draws it as a narrower line over the aligned on the diagonal). AlphaFold DB gives the overlap a third color;
  Seq\* channels have two tones, so the marks do not (a per-channel color in Seq\* would lift this).
- **The rest greyed** (`dimRest`): the structure's other residues in a neutral grey (`MarkRest`, per theme) under the
  marks, as AlphaFold DB greys what its PAE box does not cover; the marks then stand out of any coloring (a chain
  colored in the accent's hue, a pLDDT theme).
- **View state**: never a history step, never in the state tree, a snapshot or a session. A reload drops them.

## 2. Setting them

`highlight.mark({ owner, marks: [{ id, tone, label?, target }], structure?, dimRest? })` resolves each target
([targets.md](targets.md)) and replaces the owner's marks; a target that selects nothing fails (`not-found`), as
`highlight.set` does. `highlight.unmark({ owner })` clears them. Both are `ui-only` (a view's own state on the scene;
the agent highlights or selects) and take no history step.

Beside them, **the hover**: `highlight.set({ target, sequence: true })` lights the residues in the 3D view as any
highlight and in the sequence as the 3D view's hover is (`marks.hover`). A plot's cell under the pointer is that.

## 3. Showing them

- **3D**: an overpaint per representation of the marked structures, the marks' layers over the representation's own
  overpaint (its style layers'), the grey first when dimmed. Put back after Mol\* re-applies a representation's state
  (its overpaint cells changed, a representation was made), once per frame at most; the representation's own overpaint
  comes back when the marks go.
- **Sequence**: a highlight channel per mark, `marks:<owner>:<id>`, in its tone's style, on the residues it has in the
  structure the sequence shows; the hover through Seq\*'s `hover.setLoci(…, '3d')`, so it is not sent back to the 3D
  view, which lit it.
- **The view itself** draws its marks in `marks.colors` (the PAE plot: the box, its ranges on the diagonal, the edges
  projected dashed to the diagonal, the legend's swatches).
