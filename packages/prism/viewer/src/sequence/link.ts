/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import { StructureLink } from '@molstar/seq-molstar-poc';
import { linkViews, type LinkViewsOptions, type StructureView3D } from '@molstar/seq-molstar-poc/viewer';
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
  options?: LinkViewsOptions,
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
    };
  });
  return () => {
    sub.unsubscribe();
    drop();
  };
}
