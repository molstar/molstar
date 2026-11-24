/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { familySources } from './sources.js';

// The Globins for the query stories (Features/Query): the whole Pfam seed read from its file (`/data/globins`, a row form
// names its rows by their labels once they are held: it is shown first) over four of its rows inline (`/data/four`,
// named at once).

/** Four rows of the seed, as its file has them: an α and a β chain, a myoglobin and a leghemoglobin. */
export const FourRows = {
  label: ['HBA1_BOSGR/26-136', 'HBB0_MOUSE/26-142', 'MYG_CYPCA/23-136', 'LGB1_LUPLU/28-146'],
  letters: [
    'AEALERMFLSFPTTKTYFPHF......DLSQGSAQVKGHGAKVA.AA...LTKAVEHLDDLP......GALSELSDLHA..HKLRVDP....VNFKLL.SHSLLVTLASHLP......SDFTPAVHASLDKFLA.NVSTVL',
    'GETLGRLLIVYPWTQRFFDKFGNLSSAQAIMGNPRIKAHGKKVL.TS...LGLAVKNMDNLK......ETFAHLSELHC..DKLHADP....ENFKLL.GNMLVIVLSSYFG......KEFTAEAQAAWQKLVV.GVATAL',
    'GEVLTRLFKQHPETQKLFPKFVGI.ASNELAGNAAVKAHGATVL.KK...LGELLKARGDHA......AILKPLATTHA..NTHKIAL....NNFRLI.TEVLVKVMAEKAG........LDAGGQSALRRVMD.VVIGDI',
    'HRFFTLVLEIAPGAKDLFSFLKG..SSEVPQNNPDLQAHAGKVF.KL...TYEAAIQLQVNGAVA.SDATLKSLGSVHV...SKGVVD....AHFPVV.KEAILKTIKEVVG......DKWSEELNTAWTIAYD.ELAIII',
  ],
} as const;

/** The seed from its file and the four rows inline, each a section of its rows. */
export function build(): Node {
  const b = createSeqSpec({ title: 'Globins: rows and columns addressed' });
  const four = b.data.alignment('four', {
    alphabet: 'protein',
    rows: { label: [...FourRows.label], letters: [...FourRows.letters] },
    reference: FourRows.label[0],
  });
  const { aln } = familySources(b);
  b.views.view('main', (v) => {
    const f = v.section('family', { of: aln });
    f.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    f.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln });
    const s = v.section('four', { of: four });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: four });
  });
  return b.build();
}
