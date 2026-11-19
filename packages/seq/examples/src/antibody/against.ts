/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { Data, theDesign } from './data.js';
import { Chains } from './schemes.js';

// The design against the wild type (EC-29): each chain's design and wild type as the two rows of an alignment, column
// by column by their IMGT positions, so CDR-H3's two inserted residues are columns the wild type has a gap at and
// CDR-L1's deleted residue one the design has a gap at; a hovered column lights in both rows. The language has no map
// with indels between two sequences yet (`segments`, `pairwise`), so a caret on the wild type at an insertion and a
// design's CDR selected onto the wild type wait for it.

export const title = 'Antibody: the design against the wild type';

/** The rows' labels: the design's, the wild type's. */
export const RowLabels = ['design', 'trastuzumab'] as const;

export function build(d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const des = theDesign(d);
  b.views.view('main', (v) => {
    for (const c of Chains) {
      const aln = b.data.alignment(`${c}-pair`, {
        alphabet: 'protein',
        rows: { label: [...RowLabels], letters: [des[c].aligned.design, des[c].aligned.wildType] },
        reference: RowLabels[1],
      });
      const s = v.section(c, { of: aln });
      s.track('ruler', { label: c === 'light' ? 'Light chain' : 'Heavy chain' }).layer('ruler', { type: 'ruler' });
      // two rows a band: the section no taller than they are (`fill`, the default, takes the view's room left)
      s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln, height: 2 });
    }
  });
  return b.build();
}
