/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { Extents } from '@molstar/seq-core/model';
import { StructureLink } from '@molstar/seq-molstar-poc';
import { coordinatesOf, linkViews, type LinkViewsOptions, type StructureView3D } from '@molstar/seq-molstar-poc/viewer';
import type { MarksSource } from '../sources/marks.js';
import { lociSize } from '@molstar/seq-viewer';
import type { SequenceViewer } from './model.js';

/** The 3D view the sequence links to: Seq*'s contract and, optionally, the camera's reset to the whole structure. */
export interface SequenceView3D extends StructureView3D {
  /** The default framing of the structure (below the sequence band), animated; absent: the camera is left alone. */
  resetCamera?(): void;
}

/** Selection origins of the person's own moves in the sequence view (the rest: 3D, data, API, replays). */
const FromSequence: ReadonlySet<string> = new Set(['pointer', 'keyboard']);

/**
 * The two-way link between a viewer and the 3D view of the structure it shows (prism-sequence.md §2.3): hover and
 * selection both ways, focus. Made when the viewer has a structure and a 3D view of it, dropped and made again when
 * the structure changes. Returns the unlink: the control calls it when it goes, so nothing stays subscribed. The replay
 * rule makes remounting safe (linking sends the current selection as a selection, never a clear). A selection the person
 * empties in the sequence view (a click or Shift+drag removing the last stretch, Escape) also resets the 3D camera, once
 * (prism-sequence.md §9); an emptying from 3D, the selection bar, the API or a data edit leaves the camera.
 */
export function linkSequence(
  viewer: SequenceViewer,
  view3d: (structure: Structure) => SequenceView3D | undefined,
  options?: LinkViewsOptions & { readonly marks?: MarksSource },
): () => void {
  const seq = viewer.seq;
  let unlink: (() => void) | undefined;
  let link: StructureLink | undefined;
  const drop = () => {
    unlink?.();
    link?.detach();
    unlink = link = undefined;
  };
  const sub = viewer.structure.subscribe((shown) => {
    drop();
    const set = seq.graph.assembled.value?.set;
    const view = shown && view3d(shown.structure);
    if (!shown || !set || !view) return;
    link = StructureLink.attach(seq.systems, '3d', shown.structure, set, shown.complex);
    const stop = linkViews(seq, view, link, options);
    const stopMarks = options?.marks ? showMarks(seq, link, shown.structure, options.marks) : undefined;
    // Whether the selection held something: the reset follows only the move that empties it, not an empty selection again.
    let held = false;
    const sub = seq.selection.current.subscribe((sel) => {
      const parts = [sel.seq, sel.aligned, sel.columns, sel.rows, sel.items, sel.boundaries] as Parameters<
        typeof lociSize
      >[0][];
      const empty = parts.every((l) => {
        const size = lociSize(l);
        return !size || (size.items === 0 && size.length === 0);
      });
      const emptied = held && empty;
      held = !empty;
      if (emptied && FromSequence.has(sel.origin)) view.resetCamera?.();
    });
    unlink = () => {
      sub.unsubscribe();
      stop();
      stopMarks?.();
    };
  });
  return () => {
    sub.unsubscribe();
    drop();
  };
}

/**
 * A view's marks and hover in the sequence (marks.md): each mark a highlight channel (`marks:<owner>:<id>`) in its
 * tone's style (Seq*'s `highlight` or `accent`, the skin's colors, as the 3D view paints them), on the residues it has in
 * the structure shown; the hover as the 3D view's hover is (so it is not sent back to the 3D view, which lit it).
 * Returns the stop, which clears what it set.
 */
function showMarks(seq: SequenceViewer['seq'], link: StructureLink, structure: Structure, marks: MarksSource) {
  const scratch = Extents.create(16);
  const coordsOf = (loci: StructureElement.Loci) => {
    const here = StructureElement.Loci.remap(loci, structure);
    if (StructureElement.Loci.isEmpty(here)) return undefined;
    scratch.count = 0;
    return coordinatesOf(link, link.fromLoci(here, scratch));
  };
  let channels = new Set<string>();
  const marked = marks.marks.subscribe((list) => {
    const next = new Set<string>();
    for (const m of list) {
      const coords = coordsOf(m.loci);
      if (!coords) continue;
      const channel = `marks:${m.owner}:${m.id}`;
      seq.highlight.set(channel, coords, { style: m.tone });
      next.add(channel);
    }
    for (const c of channels) if (!next.has(c)) seq.highlight.clear(c);
    channels = next;
  });
  let hovering = false;
  const hovered = marks.hover.subscribe((h) => {
    const coords = h ? coordsOf(h.loci) : undefined;
    if (!coords && !hovering) return;
    hovering = !!coords;
    seq.hover.setLoci(coords ?? null, '3d');
  });
  return () => {
    marked.unsubscribe();
    hovered.unsubscribe();
    for (const c of channels) seq.highlight.clear(c);
    if (hovering) seq.hover.setLoci(null, '3d');
  };
}
