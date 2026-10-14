/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fixture of the text-only letters row's marks (`letters-marks.test.ts`, `letters-marks.heap.ts`).

import { createSeqViewSpec } from '@molstar/seq-core/builder';

const UNIT = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV';
/** 400 residues: the whole row is a device px or two a position at the widest fit. */
const WT = UNIT.repeat(5);
const NO_CELL_COLOR = { theme: 'rules', params: { rules: [] } } as const;

/** The wild type with a design under it missing residues 30–39, its own letters with or without a cell color. */
export function marksDoc(color: boolean) {
  const b = createSeqViewSpec();
  const wt = b.data.sequence({ ref: 'wt', alphabet: 'protein', value: WT });
  const design = b.data.sequence({ ref: 'd', alphabet: 'protein', value: WT.slice(0, 30) + WT.slice(40) });
  b.data.mapping({
    ref: 'd.on_wt',
    from: wt,
    to: design,
    method: 'pairwise',
    gapped: [WT, `${WT.slice(0, 30)}${'-'.repeat(10)}${WT.slice(40)}`],
  });
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ sequence: wt });
  v.track({ label: 'Design' }).letters({ of: design, ...(color ? {} : { color: NO_CELL_COLOR }) });
  v.track({ label: 'Plain' }).letters(color ? {} : { color: NO_CELL_COLOR });
  return b.getDocument();
}
