/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, modified, type SeqSpecBuilder } from '@molstar/seqx-core/lang';

// Molecules drawn (the design's Molecules: On an axis, A section of a complex): a small generated glycoprotein complex,
// written inline: an Fc-like chain in copies with modified residues (MSE on M1, SEP on S21, a three-letter CRO over
// T31 Y32 G33), a glycan on N47 of the first chain (the sequon N-L-S), a disulfide from C5 of the first chain to C71 of
// the second, ATP with its heavy atoms, a sodium ion bound to D64's OD1 and ATP's γ-phosphate, a water (left out by
// default), and a tRNA-like sequence with modified nucleotides (2MG, PSU, H2U | U). Not an entry: the shapes the
// design names, at a size a story reads at a glance (the recorded entries are the examples' DT, DC and DG).

const FcLetters = 'MKTACIAKQRQISFVKSHFSSQLEERLGLITYGAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVCSLAKWKRQTLGQ';
const TrnaLetters = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA';

interface Molecules {
  readonly fc: string;
  readonly trna: string;
  readonly complex: string;
  readonly links: string;
  /** The units by id: the chains `A`, `B`, …, then `glycan`, `atp`, `na`, `water`. */
  readonly units: Readonly<Record<string, string>>;
  readonly chains: readonly string[];
}

/** The chains' unit ids: `A`, `B`, … then `C26`, `C27`, … */
export function chainIds(n: number): string[] {
  return Array.from({ length: n }, (_, i) => (i < 26 ? String.fromCharCode(65 + i) : `C${i}`));
}

/** The dataset written into a builder: `chains` copies of the Fc chain (2 by default). */
export function moleculesData(b: SeqSpecBuilder, chains = 2): Molecules {
  const fc = b.data.sequence('fc', {
    alphabet: 'protein',
    name: 'Fc',
    letters: FcLetters,
    modified: modified([
      { position: 0, component: 'MSE' },
      { position: 20, component: 'SEP' },
      { position: 30, component: 'CRO', length: 3 },
    ]),
  });
  const glycan = b.data.molecule('glycan', {
    class: 'branched',
    name: 'Fc glycan',
    residues: ['NAG', 'NAG', 'BMA', 'MAN', 'MAN'],
    identifiers: identifiers({ glytoucan: 'G15407YE' }),
  });
  const atp = b.data.molecule('atp', {
    residues: ['ATP'],
    atoms: {
      residue: [0, 0, 0, 0, 0, 0, 0, 0],
      name: ['PG', 'O1G', 'O2G', 'PB', 'PA', 'N1', 'C2', 'N3'],
      element: ['P', 'O', 'O', 'P', 'P', 'N', 'C', 'N'],
    },
    identifiers: identifiers([
      ['ccd', 'ATP'],
      ['smiles', 'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O'],
    ]),
  });
  const na = b.data.molecule('na', { class: 'ion', residues: ['NA'], identifiers: identifiers([['ccd', 'NA']]) });
  const water = b.data.molecule('water', { class: 'water', residues: ['HOH'] });
  b.data.compound('SEP', { name: 'phosphoserine', parents: ['SER'] });
  b.data.compound('MSE', { name: 'selenomethionine', parents: ['MET'] });
  b.data.compound('CRO', { name: 'chromophore', parents: ['THR', 'TYR', 'GLY'] });
  b.data.compound('PSU', { name: 'pseudouridine-5′-monophosphate', parents: ['U'] });
  b.data.compound('2MG', { name: '2N-methylguanosine-5′-monophosphate', parents: ['G'] });
  const trna = b.data.sequence('trna', {
    alphabet: 'rna',
    name: 'tRNA',
    letters: TrnaLetters,
    modified: modified([
      { position: 9, component: '2MG' },
      { position: 16, component: 'H2U' },
      { position: 16, component: 'U' },
      { position: 38, component: 'PSU' },
      { position: 54, component: 'PSU' },
    ]),
  });
  const ids = chainIds(chains);
  const units: Record<string, string> = {};
  let links = '';
  const complex = b.data.complex('complex', { name: 'Fc with its glycan' }, (c) => {
    for (const id of ids) units[id] = c.unit(id, { of: fc, chain: id, label: `Fc ${id}` });
    units.glycan = c.unit('glycan', { of: glycan, chain: ids.length < 26 ? 'X' : 'GLC' });
    units.atp = c.unit('atp', { of: atp, chain: 'Y' });
    units.na = c.unit('na', { of: na, chain: 'Z' });
    units.water = c.unit('water', { of: water, chain: 'W' });
    links = c.links('links', {
      rows: {
        a_unit: ['A', 'A', 'na', 'na'],
        a_start: [46, 4, 0, 0],
        a_end: [47, 5, 1, 1],
        b_unit: ['glycan', ids[1] ?? 'A', 'A', 'atp'],
        b_start: [0, 70, 63, 0],
        b_end: [1, 71, 64, 1],
        type: ['covalent', 'disulfide', 'metal', 'metal'],
        a_atom: ['ND2', 'SG', 'NA', 'NA'],
        b_atom: ['C1', 'SG', 'OD1', 'O1G'],
      },
    });
  });
  return { fc, trna, complex, links, units, chains: ids };
}

/** The complex as a section of it: a block a unit, the chains merged `×2`, the molecules as chips; links as arcs. */
export function complexView(): Node {
  const b = createSeqSpec({ title: 'Molecules: a section of a complex' });
  const d = moleculesData(b);
  b.views.view('main', (v) =>
    v.section('complex', { of: d.complex as never }, (s) => {
      s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
      s.track('links', { label: 'Links' }).layer('links', { type: 'pairs', data: d.links as never });
    }),
  );
  return b.build();
}

/** The forms side by side: the molecules whole, then the glycan's residues and ATP's atoms (explicit blocks). */
export function formsView(): Node {
  const b = createSeqSpec({ title: 'Molecules: whole, residues, atoms' });
  const d = moleculesData(b);
  b.views.view('main', (v) => {
    v.section('whole', { of: d.complex as never, classes: ['branched', 'ligand', 'ion'] }, (s) =>
      s.track('letters', { label: 'Whole' }).layer('letters', { type: 'letters' }),
    );
    v.section('parts', {}, (s) => {
      s.block('glycan', { of: d.units.glycan as never, form: 'residues', label: 'Fc glycan: residues' });
      s.block('atp', { of: d.units.atp as never, form: 'atoms', label: 'ATP: atoms' });
      s.track('letters', { label: 'Residues, atoms' }).layer('letters', { type: 'letters' });
    });
  });
  return b.build();
}

/** A modified residue drawn four ways: `mark` (the default), `code`, `tint`, `none`. */
export function modifiedView(): Node {
  const b = createSeqSpec({ title: 'Molecules: modified residues' });
  const d = moleculesData(b);
  b.views.view('main', (v) => {
    for (const m of ['mark', 'code', 'tint', 'none'] as const)
      v.section(m, { of: d.trna as never }, (s) =>
        s.track('letters', { label: `modified: ${m}` }).layer('letters', { type: 'letters', modified: m }),
      );
    v.section('fc', { of: d.fc as never }, (s) =>
      s.track('letters', { label: 'Fc (MSE, SEP, CRO)' }).layer('letters', { type: 'letters' }),
    );
  });
  return b.build();
}

/** The complex as text: the chains' letters, then the molecules' chips flowing like words (EC-99). */
export function textView(): Node {
  const b = createSeqSpec({ title: 'Molecules: as text' });
  const d = moleculesData(b);
  b.views.view('main', { layout: 'text' }, (v) =>
    v.section('complex', { of: d.complex as never, copies: 'each' }, (s) =>
      s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' }),
    ),
  );
  return b.build();
}

/** A capsid-like complex: 200 copies of the chain merged into one block `×200`, its glycan a chip (EC-67). */
export function copiesView(): Node {
  const b = createSeqSpec({ title: 'Molecules: 200 copies' });
  const d = moleculesData(b, 200);
  b.views.view('main', (v) =>
    v.section('complex', { of: d.complex as never }, (s) => {
      s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
      s.track('links', { label: 'Links' }).layer('links', { type: 'pairs', data: d.links as never });
    }),
  );
  return b.build();
}
