/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Boundaries over their track's letters, for the pins' test, conformance and heap guard: case 15's design row (the
// letters of a design through its mapping, its insertions' marks over them) and a block's own letters with the same
// marks, at the block's start, its middle (20) and its end.

import { createSeqViewSpec } from '@molstar/seq-core/builder';

export const S = 'v.main';
export const WT = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV';
export const N = WT.length;
/** The design: WT with GSGS after 20, as a design row's letters (`of`). */
const D = `${WT.slice(0, 20)}GSGS${WT.slice(20)}`;
const MARK = '#123456';

export function pinsDoc(height = 1) {
  const b = createSeqViewSpec();
  const wt = b.data.sequence({ ref: 'wt', value: WT });
  const d = b.data.sequence({ ref: 'd', value: D });
  b.data.mapping({
    ref: 'd.on_wt',
    from: wt,
    to: d,
    method: 'pairwise',
    gapped: [`${WT.slice(0, 20)}----${WT.slice(20)}`, D],
  });
  const marks = b.data.annotation({
    ref: 'ins',
    geometry: 'ranges',
    target: wt,
    value: { boundary: [0, 20, N], label: ['N-terminal', 'graft', 'C-terminal'] },
  });
  const sec = b.view({ ref: 'v', label: 'V' }).section({ ref: S }).block({ sequence: wt });
  const ranges = {
    data: marks,
    color: { theme: 'uniform', params: { color: MARK } },
    glyph: { theme: 'uniform', params: { glyph: { shape: 'triangle', label: 'insertion' } } },
    radius: { theme: 'uniform', params: { size: 9 } },
  } as never;
  // A design row (case 15's: letters of the design, its marks over them), and the block's own letters.
  sec.track({ ref: 't-of', label: 'Design', height }).letters({ of: d }).ranges(ranges);
  sec.track({ ref: 't-own', label: 'WT', height }).letters().ranges(ranges);
  return b.getDocument();
}
