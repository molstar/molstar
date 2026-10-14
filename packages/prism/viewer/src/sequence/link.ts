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
import type { SequenceViewer } from './model.js';

/**
 * The two-way link between a viewer and the 3D view of the structure it shows (prism-sequence.md §2.3): hover and
 * selection both ways, focus. Made when the viewer has a structure and a 3D view of it, dropped and made again when
 * the structure changes. Returns the unlink: the control calls it when it goes, so nothing stays subscribed. The replay
 * rule makes remounting safe (linking sends the current selection as a selection, never a clear).
 */
export function linkSequence(
  viewer: SequenceViewer,
  view3d: (structure: Structure) => StructureView3D | undefined,
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
    unlink = linkViews(seq, view, link, options);
  });
  return () => {
    sub.unsubscribe();
    drop();
  };
}
