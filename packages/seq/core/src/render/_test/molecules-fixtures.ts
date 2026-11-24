/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A complex drawn (_spec/README.md "Render", Blocks; Layer types): an Fc-like chain in copies, a glycan on the first,
// ATP with its atoms, a water and a sodium ion, the links between them; a tRNA-like sequence with modified nucleotides.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, modified } from '@molstar/seqx-core/lang';

export const FcLetters = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ';

export interface MoleculesOptions {
  readonly copies?: 'merge' | 'each';
  readonly classes?: readonly string[];
  /** Copies of the Fc chain (2 by default). */
  readonly chains?: number;
  /** The letters layer's `modified`. */
  readonly modified?: 'mark' | 'code' | 'tint' | 'none';
  readonly layout?: 'tracks' | 'text';
  /** A text's numbers (`numbers`). */
  readonly numbers?: 'start' | 'end' | 'both';
  /** ATP's block form (by an explicit block list: the complex's units in order, ATP's form as given). */
  readonly atpForm?: 'whole' | 'residues' | 'atoms';
  /** Overrides of derived units: `block` children by a unit's id, with no `of`. */
  readonly overrides?: Readonly<Record<string, { form?: 'whole' | 'residues' | 'atoms'; label?: string }>>;
  /** A second sodium ion's unit (`n2`), after the others. */
  readonly extraIon?: boolean;
  /** The links' arcs (a pairs track; true but in text). */
  readonly links?: boolean;
  /** The links' sites: a ranges track over the links, each end a mark labeled by the other. */
  readonly sites?: boolean;
  /** Pairs within the Fc chain on its molecule (no unit named): a pairs track, drawn on each copy split. */
  readonly within?: boolean;
  /** The section's letters layer's `color` (a theme form). */
  readonly color?: Node;
  /** A second letters layer, `second`, on the same track (a track of two layers: FX29's layer names). */
  readonly second?: Node;
}

/** The Fc chains' unit ids: `A`, `B`, … */
export function chainIds(n: number): string[] {
  return Array.from({ length: n }, (_, i) => (i < 26 ? String.fromCharCode(65 + i) : `X${i}`));
}

/** A complex with polymers, a glycan, ATP, an ion and water, and a section of it; a modified tRNA beside it. */
export function moleculesSpec(o: MoleculesOptions = {}): Node {
  const b = createSeqSpec({ title: 'molecules' });
  const fc = b.data.sequence('fc', {
    alphabet: 'protein',
    letters: FcLetters,
    modified: modified([
      { position: 0, component: 'MSE' },
      { position: 20, component: 'SEP' },
      { position: 30, component: 'CRO', length: 3 },
    ]),
  });
  const glycan = b.data.molecule('glycan', {
    class: 'branched',
    residues: ['NAG', 'NAG', 'BMA', 'MAN', 'MAN'],
    identifiers: identifiers({ glytoucan: 'G15407YE' }),
  });
  const atp = b.data.molecule('atp', {
    residues: ['ATP'],
    atoms: { residue: [0, 0, 0, 0], name: ['PG', 'O1G', 'N1', 'C2'], element: ['P', 'O', 'N', 'C'] },
    identifiers: identifiers([['ccd', 'ATP']]),
  });
  const na = b.data.molecule('na', { class: 'ion', residues: ['NA'] });
  const hoh = b.data.molecule('hoh', { class: 'water', residues: ['HOH'] });
  b.data.compound('SEP', { name: 'phosphoserine', parents: ['SER'] });
  const trna = b.data.sequence('trna', {
    alphabet: 'rna',
    letters: 'GCGGAUUUAGCUCAGUUGGGAGAGC',
    modified: modified([
      { position: 9, component: '2MG' },
      { position: 15, component: 'H2U' },
      { position: 15, component: 'U' },
    ]),
  });
  const within = o.within
    ? b.data.pairs('within', {
        a: fc,
        b: fc,
        rows: { a_start: [4, 10, 20], a_end: [5, 11, 21], b_start: [70, 40, 30], b_end: [71, 41, 31] },
      })
    : undefined;
  b.data.compound('2MG', { name: '2N-methylguanosine-5’-monophosphate', parents: ['G'] });
  const ids = chainIds(o.chains ?? 2);
  const units: Record<string, string> = {};
  b.data.complex('c', { name: 'Fc with a glycan' }, (c) => {
    for (const id of ids) units[id] = c.unit(id, { of: fc, chain: id });
    units.w = c.unit('w', { of: hoh, chain: 'W' });
    units.g = c.unit('g', { of: glycan, chain: 'G', label: 'Glycan' });
    units.l = c.unit('l', { of: atp, chain: 'L' });
    units.n = c.unit('n', { of: na, chain: 'N' });
    if (o.extraIon) units.n2 = c.unit('n2', { of: na, chain: 'M', label: 'NA 2' });
    c.links('links', {
      rows: {
        a_unit: ['A', 'A', 'l'],
        a_start: [10, 4, 0],
        a_end: [11, 5, 1],
        b_unit: ['g', ids[1] ?? 'A', 'A'],
        b_start: [0, 4, 50],
        b_end: [1, 5, 51],
        type: ['covalent', 'disulfide', 'metal'],
      },
    });
  });
  b.views.view('v', (v) => {
    v.section(
      'complex',
      {
        of: '/data/c',
        ...(o.copies ? { copies: o.copies } : {}),
        ...(o.classes ? { classes: [...o.classes] } : {}),
        ...(o.layout ? { layout: o.layout } : {}),
        ...(o.numbers ? { numbers: o.numbers } : {}),
      },
      (s) => {
        for (const [id, p] of Object.entries(o.overrides ?? {})) s.block(id, p);
        if (o.atpForm)
          for (const id of [...ids, 'g', 'l', 'n'])
            s.block(id, { of: units[id] as never, ...(id === 'l' ? { form: o.atpForm } : {}) });
        const letters = s.track('letters');
        letters.layer('letters', {
          type: 'letters',
          ...(o.modified ? { modified: o.modified } : {}),
          ...(o.color ? { color: o.color } : {}),
        });
        if (o.second) letters.layer('second', { type: 'letters', color: o.second });
        if (o.sites) s.track('sites').layer('sites', { type: 'ranges', data: '/data/c/links' as never });
        if (within) s.track('within').layer('within', { type: 'pairs', data: within });
        if (o.layout !== 'text' && o.links !== false)
          s.track('links').layer('links', { type: 'pairs', data: '/data/c/links' as never });
      },
    );
    v.section('trna', { of: trna }, (s) =>
      s.track('t').layer('letters', { type: 'letters', ...(o.modified ? { modified: o.modified } : {}) }),
    );
  });
  return b.build();
}
