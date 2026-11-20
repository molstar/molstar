/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/seqx-core/query';
import { nearestOffsetOf, type Projector } from '@molstar/seqx-viewer/models';
import type { SeqViewer } from './viewer.js';

// What a selection moves beyond its own pass (_spec/README.md "Composition"): the sections that show none of a
// selection made elsewhere scrolled to it (`reveal.linked`), and its dynamic nodes held in the store while it is
// selected (`selection.holdNodes`).

/**
 * A selection made elsewhere (another section, a link, a structure) scrolls each section that shows none of it to its
 * first part (_spec/README.md "Composition", `reveal.linked`): text sections by default, as a move by hand (no
 * step of the history). Not the section the selection is being made in (a drag's, the cursor's): it shows it already.
 */
export function revealLinked(viewer: SeqViewer, projector: Projector): void {
  const mode = viewer.config.reveal.linked;
  const resolved = viewer.view.resolved.value;
  const shown = viewer.selection.value.value.shown;
  if (mode === 'none' || !resolved || viewer.retired || !shown.parts.length) return;
  const own = viewer.input.selectingIn ?? viewer.selection.cursor.value?.section;
  let loc: Location | undefined;
  for (const s of resolved.sections) {
    if (s.address === own || (mode === 'text' && s.layout !== 'text')) continue;
    const v = viewer.view.get(s.address);
    if (!v || v.pending) continue;
    loc ??= projector.of(shown, { own: true }) ?? shown;
    const r = loc.parts.find((p) => p.system === s.axis.system)?.ranges;
    if (!r?.length) continue;
    let first = Number.POSITIVE_INFINITY;
    let shownNow = false;
    for (let i = 0; i + 1 < r.length && !shownNow; i += 2) {
      const a = nearestOffsetOf(s.axis, r[i]);
      const b = nearestOffsetOf(s.axis, r[i + 1]);
      if (a === undefined || b === undefined) continue;
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      if (hi > v.start && lo < v.end) shownNow = true;
      first = Math.min(first, lo);
    }
    if (!shownNow && Number.isFinite(first)) viewer.view.reveal(s.address, first, first + 1);
  }
}

/** The selection's dynamic nodes held while it is selected (`selection.holdNodes`). */
export function holdSelection(viewer: SeqViewer): void {
  if (!viewer.config.selection.holdNodes || viewer.retired) return;
  const sel = viewer.selection.value.value;
  const addresses: string[] = [];
  for (const k of sel.elements.keys) addresses.push(k.owner);
  for (const part of sel.shown.parts) addresses.push(part.system);
  viewer.store.holdSession('selection', viewer.spec.model.value, addresses);
}
