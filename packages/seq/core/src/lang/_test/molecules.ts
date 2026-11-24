/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures of the molecule kinds (_spec/README.md, Language, Molecules and complexes): a glycoprotein with a ligand,
// written with every new param, and one misuse a `seq.lang.molecule.*`, `modified.*`, `unit.*`, `link.*` code.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, modified } from '../builder.js';
import type { LangIssues } from '../issues.js';

/** 1EHZ-like letters with modified nucleotides, an Fc chain with its N297-like glycan, ATP with atoms. */
export function moleculesSpec(): Node {
  const b = createSeqSpec({ title: 'molecules' });
  const host = b.data.source('host', { format: 'mmcif', uri: 'https://files.example/1abc.cif', ref: 'model:7' });
  const trna = b.data.sequence('trna', {
    alphabet: 'rna',
    letters: 'GCGGAUUUAGCUCAG',
    modified: modified([
      { position: 9, component: '2MG' },
      { position: 12, component: 'H2U' },
      { position: 12, component: 'U' },
    ]),
  });
  const fc = b.data.sequence('fc', {
    alphabet: 'protein',
    letters: 'MKTAYIAKQRQISFVKSHFSRQ',
    modified: modified([
      { position: 0, component: 'MSE' },
      { position: 3, component: 'CRO', length: 3 },
    ]),
  });
  const glycan = b.data.molecule(
    'glycan',
    {
      class: 'branched',
      residues: ['NAG', 'NAG', 'BMA'],
      identifiers: identifiers({ glytoucan: 'G15407YE', wurcs: 'WURCS=2.0/3,3,2/…' }),
    },
    (m) =>
      m.links('bonds', {
        rows: { a_start: [1, 2], a_end: [2, 3], b_start: [0, 1], b_end: [1, 2], type: ['glycosidic', 'glycosidic'] },
      }),
  );
  const atp = b.data.molecule('atp', {
    name: 'ATP',
    residues: ['ATP'],
    atoms: { residue: [0, 0, 0], name: ['PG', 'O1G', 'N1'], element: ['P', 'O', 'N'] },
    identifiers: identifiers([
      ['ccd', 'ATP'],
      ['smiles', 'c1nc(c2c(n1)n(cn2)…)N'],
    ]),
  });
  b.data.compound('CRO', { name: 'chromophore', parents: ['THR', 'TYR', 'GLY'], class: 'L-PEPTIDE LINKING' });
  b.data.compound('2MG', { name: '2N-methylguanosine-5’-monophosphate', parents: ['G'] });
  let unitA = '';
  b.data.complex('c', { name: 'Fc with a glycan' }, (c) => {
    unitA = c.unit(
      'A',
      { of: fc, range: [0, 20], observed: { start: [2], end: [18] }, chain: 'A', auth_chain: 'H', operator: '1_555' },
      (u) =>
        u.numbering('auth', {
          method: 'labels',
          scheme: 'auth',
          labels: Array.from({ length: 20 }, (_, i) => (i < 2 ? null : `${i + 230}`)),
        }),
    );
    c.unit('B', { of: fc, range: [0, 20], chain: 'B' });
    c.unit('G', { of: glycan, chain: 'C' });
    c.unit('L', { of: atp, chain: 'D', group: 'PRD_000001' });
    c.links('links', {
      rows: {
        a_unit: ['A', 'A'],
        a_start: [10, 4],
        a_end: [11, 5],
        b_unit: ['G', 'B'],
        b_start: [0, 4],
        b_end: [1, 5],
        type: ['covalent', 'disulfide'],
        a_atom: ['ND2', 'SG'],
        b_atom: ['C1', 'SG'],
        occupancy: [0.6, 1],
      },
    });
  });
  const sites = b.data.ranges('sites', { target: unitA, rows: { start: [10], end: [11], category: ['glycan'] } });
  b.data.values('plddt', { target: fc, source: host, field: 'plddt' });
  b.views.view('v', (v) => {
    v.section('complex', { of: '/data/c', classes: ['polymer', 'branched', 'ligand'], copies: 'each' }, (s) =>
      s.track('letters').layer('letters', { type: 'letters', modified: 'code' }),
    );
    v.section('blocks', {}, (s) => {
      s.block('fc', { of: unitA, form: 'residues' });
      s.block('atp', { of: atp, form: 'atoms' });
      s.block('glycan', { of: glycan });
      s.track('sites').layer('sites', { type: 'ranges', data: '/data/c/links' });
    });
    v.section('trna', { of: trna }, (s) => s.track('t').layer('marks', { type: 'letters', modified: 'mark' }));
    v.section('sites', { of: unitA }, (s) => s.track('t').layer('s', { type: 'ranges', data: sites }));
  });
  return b.build();
}

type MoleculeCode = Extract<
  keyof typeof LangIssues.codes,
  `seq.lang.${'molecule' | 'modified' | 'unit' | 'link'}.${string}` | 'seq.lang.source.host_only'
>;

const fcLetters = 'MKTAYIAK';

/** One spec a code: valid for the description, wrong in the language's own way. */
export const MoleculeMisuses: Record<MoleculeCode, () => Node> = {
  'seq.lang.source.host_only': () => {
    const b = createSeqSpec();
    b.data.source('s', { format: 'mmcif', ref: 'edited-model' });
    return b.build();
  },
  'seq.lang.modified.order': () => {
    const b = createSeqSpec();
    b.data.sequence('s', {
      alphabet: 'protein',
      letters: fcLetters,
      modified: modified([
        { position: 4, component: 'MSE' },
        { position: 1, component: 'SEP' },
      ]),
    });
    return b.build();
  },
  'seq.lang.modified.bounds': () => {
    const b = createSeqSpec();
    b.data.sequence('s', {
      alphabet: 'protein',
      letters: fcLetters,
      modified: modified([{ position: 6, component: 'CRO', length: 3 }]),
    });
    return b.build();
  },
  'seq.lang.modified.length': () => {
    const b = createSeqSpec();
    b.data.sequence('s', {
      alphabet: 'protein',
      letters: fcLetters,
      modified: modified([
        { position: 1, component: 'CRO', length: 3 },
        { position: 2, component: 'SEP' },
      ]),
    });
    return b.build();
  },
  'seq.lang.molecule.content': () => {
    const b = createSeqSpec();
    b.data.molecule('m', { name: 'nothing' });
    return b.build();
  },
  'seq.lang.molecule.atoms': () => {
    const b = createSeqSpec();
    b.data.molecule('m', { residues: ['ATP'], atoms: { residue: [0, 1], name: ['PG', 'X'], element: ['P', 'O'] } });
    return b.build();
  },
  'seq.lang.molecule.hydrogen': () => {
    const b = createSeqSpec();
    b.data.molecule('m', { residues: ['HOH'], atoms: { residue: [0, 0], name: ['O', 'H1'], element: ['O', 'H'] } });
    return b.build();
  },
  'seq.lang.molecule.identifier': () => {
    const b = createSeqSpec();
    b.data.compound('ATP', { identifiers: identifiers({ 'CCD:x': 'ATP' }) });
    return b.build();
  },
  'seq.lang.molecule.form': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: fcLetters });
    b.views.view('v', (v) => v.section('x', {}, (sec) => sec.block('b', { of: s, form: 'atoms' })));
    return b.build();
  },
  'seq.lang.unit.range': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: fcLetters });
    b.data.complex('c', {}, (c) => c.unit('A', { of: s, range: [6, 2] }));
    return b.build();
  },
  'seq.lang.unit.observed': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: fcLetters });
    b.data.complex('c', {}, (c) => c.unit('A', { of: s, range: [2, 6], observed: { start: [1], end: [5] } }));
    return b.build();
  },
  'seq.lang.link.target': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: fcLetters });
    b.data.complex('c', {}, (c) => {
      c.unit('A', { of: s });
      c.links('l', {
        a: s,
        rows: { a_unit: ['A'], a_start: [0], a_end: [1], b_unit: ['A'], b_start: [4], b_end: [5] },
      } as never);
    });
    return b.build();
  },
  'seq.lang.link.unit': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: fcLetters });
    b.data.complex('c', {}, (c) => {
      c.unit('A', { of: s });
      c.links('l', { rows: { a_unit: ['A'], a_start: [0], a_end: [1], b_unit: ['Z'], b_start: [4], b_end: [5] } });
    });
    return b.build();
  },
  'seq.lang.link.bounds': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: fcLetters });
    b.data.complex('c', {}, (c) => {
      c.unit('A', { of: s, range: [0, 4] });
      c.links('l', { rows: { a_unit: ['A'], a_start: [0], a_end: [1], b_unit: ['A'], b_start: [5], b_end: [6] } });
    });
    return b.build();
  },
  'seq.lang.link.position': () => {
    const b = createSeqSpec();
    b.data.molecule('m', { residues: ['NAG', 'NAG', 'NAG'] }, (m) =>
      m.links('l', { rows: { a_start: [0], a_end: [2], b_start: [2], b_end: [3] } }),
    );
    return b.build();
  },
  'seq.lang.link.occupancy': () => {
    const b = createSeqSpec();
    b.data.molecule('m', { residues: ['NAG', 'NAG'] }, (m) =>
      m.links('l', { rows: { a_start: [1], a_end: [2], b_start: [0], b_end: [1], occupancy: [1.5] } }),
    );
    return b.build();
  },
};
