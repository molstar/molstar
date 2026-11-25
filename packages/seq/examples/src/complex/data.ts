/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { printAddress } from '@molstar/core/repr/tree';
import { DefaultOpenLimits } from '@molstar/seqx-core/data';
import { readMatrix } from '@molstar/seqx-core/data/records';
import type { SeqSpecBuilder } from '@molstar/seqx-core/lang';
import { uniprotTransform } from '../files/uniprot.js';
import type { Derived } from './derived.js';
import { DerivedData } from './derived-data.js';

// The Complex dataset: CDK2 with cyclin A. The crystal, PDB 1FIN (two copies of the complex, an ATP on each CDK2), its
// entities' letters and secondary structure read by location from `1FIN.seq`; the contacts across each copy and each
// ligand's contacting residues derived from the coordinates at recording and held here (`derived-data.ts`, a `pairs`
// node is inline rows only), compared with `contacts.seq` by the fixtures test. The ligands besides 1FIN's ATP come
// from six CC0 entries of CDK2 (AlphaFill's transplants are not taken: no CC BY 4.0 or looser licence is stated). The
// prediction: UniProt P24941 and P20248 (their letters and features) with AlphaFold DB's per-residue confidence of
// each chain, read by location. The PAE over the complex is synthetic (no permissive prediction of the complex): made
// here from a seed and 1FIN's interface (`complexPaeDocument`), AlphaFold-shaped, read through the matrix read.

export const Files = {
  molecule: 'fixtures/complex/1FIN.seq',
  contacts: 'fixtures/complex/contacts.seq',
  cdk2: 'fixtures/complex/P24941.json',
  cyclin: 'fixtures/complex/P20248.json',
  cdk2Confidence: 'fixtures/complex/AF-P24941-F1-confidence_v6.json',
  cyclinConfidence: 'fixtures/complex/AF-P20248-F1-confidence_v6.json',
} as const;

/**
 * AlphaFill's usage policy as checked when the Complex was recorded: a policy of its own (free use, attribution of
 * AlphaFill and AlphaFold, modified files only with the parent attributed), no CC BY 4.0 or looser licence stated, so
 * its transplants are not taken (the data's decisions, 1).
 */
export const AlphaFillPolicy = { url: 'https://alphafill.eu/license', checked: '2026-11-19' } as const;

/** The systems of `1FIN.seq` (and `contacts.seq`): CDK2's entity and the cyclin A construct's. */
export const CrystalSystems = { cdk2: 'cdk2', cyclin: 'cyclin-a' } as const;

/** The chains' lengths: CDK2 (P24941, whole in 1FIN), the cyclin A construct of 1FIN, cyclin A2 (P20248). */
export const Lengths = { cdk2: 298, construct: 260, cyclin: 432 } as const;

/** The derived data the views hold (the fixture by default). */
export const Data: Derived = DerivedData;

/** The construct's offset onto P20248 (its first residue's author number less one: 172). */
export const constructOffset = (d: Derived = Data) => d.cyclinFirst - 1;

/**
 * Features named in the literature, on CDK2 (0-based, half-open; UniProt P24941's numbering), with no UniProt feature
 * of their own: the PSTAIRE (αC) helix, the activation loop from DFG to APE, the T160 phosphosite (UniProt's modified
 * residue too), and cyclin A2's cyclin box folds (on P20248). Curated here, illustrative of the literature's spans.
 */
export const Cdk2Features = {
  start: [44, 144, 159],
  end: [56, 172, 160],
  id: ['pstaire', 'aloop', 't160'],
  label: ['PSTAIRE helix', 'activation loop', 'T160'],
};
export const CyclinFeatures = {
  start: [208, 308],
  end: [306, 400],
  id: ['box-n', 'box-c'],
  label: ['cyclin box (N)', 'cyclin box (C)'],
};

/** A `pairs` node's rows from contacts: CDK2's residue against cyclin A's (on the construct, or shifted onto P20248). */
export function contactRows(options: { readonly copy?: number; readonly shift?: number } = {}, d: Derived = Data) {
  const shift = options.shift ?? 0;
  const kept = d.contacts.filter((c) => options.copy === undefined || c[0] === options.copy);
  return {
    a_start: kept.map((c) => c[1]),
    a_end: kept.map((c) => c[1] + 1),
    b_start: kept.map((c) => c[2] + shift),
    b_end: kept.map((c) => c[2] + shift + 1),
    type: kept.map(() => 'contact'),
    id: kept.map((c) => `${c[0]}:${c[1] + 1}-${c[2] + 1 + shift}`),
  };
}

/** The interface: each side's contacting residues, sorted, as ranges rows (a range a run), on its own system. */
export function interfaceRows(side: 'cdk2' | 'cyclin', d: Derived = Data) {
  const at = [...new Set(d.contacts.map((c) => (side === 'cdk2' ? c[1] : c[2])))].sort((a, b) => a - b);
  const start: number[] = [];
  const end: number[] = [];
  for (const p of at) {
    if (end.length && end[end.length - 1] === p) end[end.length - 1] = p + 1;
    else {
      start.push(p);
      end.push(p + 1);
    }
  }
  return { start, end, id: start.map((s, k) => `${side}-${s + 1}-${end[k]}`), label: start.map(() => 'interface') };
}

/** A ligand of the dataset by its component (`ATP` the first copy's). */
export const ligand = (comp: string, d: Derived = Data) => d.ligands.find((l) => l.comp === comp)!;

/** The ligands of the other entries (not 1FIN's ATPs), best resolution first. */
export const entryLigands = (d: Derived = Data) =>
  d.ligands.filter((l) => l.entry !== '1FIN').sort((a, b) => a.resolution - b.resolution);

/**
 * A ligand's contacts as `pairs` rows: from its one position (`per: 'ligand'`, a pair a residue) or from each heavy
 * atom (`per: 'atom'`, a pair an atom and residue) to CDK2's residues.
 */
export function ligandContactRows(l: Derived['ligands'][number], per: 'ligand' | 'atom' = 'ligand') {
  const pairs =
    per === 'atom'
      ? l.contacts.map(([r, a]) => [a, r] as const)
      : [...new Set(l.contacts.map(([r]) => r))].sort((x, y) => x - y).map((r) => [0, r] as const);
  return {
    a_start: pairs.map((p) => p[0]),
    a_end: pairs.map((p) => p[0] + 1),
    b_start: pairs.map((p) => p[1]),
    b_end: pairs.map((p) => p[1] + 1),
    type: pairs.map(() => 'contact'),
    id: pairs.map((p) => `${l.comp}-${per === 'atom' ? `${l.atoms[p[0]]}-` : ''}${p[1] + 1}`),
  };
}

/** 1FIN's source and its two entities' sequences, their letters and secondary structure read from the file. */
export function crystal(b: SeqSpecBuilder, uri: string = Files.molecule) {
  const source = b.data.source('fin', {
    format: 'seq',
    uri,
    attribution: 'RCSB PDB 1FIN (Jeffrey et al. 1995)',
    licence: 'CC0-1.0',
  });
  const cdk2 = b.data.sequence('cdk2', { alphabet: 'protein', source, name: CrystalSystems.cdk2 });
  const cyclin = b.data.sequence('cyclin', { alphabet: 'protein', source, name: CrystalSystems.cyclin });
  /** Each entity's secondary structure (ranges with `chain` and `category`). */
  const cdk2Secondary = b.data.ranges('cdk2-secondary', { target: cdk2, source });
  const cyclinSecondary = b.data.ranges('cyclin-secondary', { target: cyclin, source });
  return { source, cdk2, cyclin, cdk2Secondary, cyclinSecondary };
}

/**
 * AlphaFold DB's per-residue confidence file (`AF-<accession>-F1-confidence_v6.json`) as a `values` node on `target`: the
 * source reads the file as published (`residueNumber`, `confidenceScore`), the node is the values over the protein. The
 * services showcase reads the same file through it.
 */
export function afdbConfidence(
  b: SeqSpecBuilder,
  o: {
    readonly source: string;
    readonly values: string;
    readonly accession: string;
    readonly uri: string;
    readonly target: string;
  },
): string {
  const source = b.data.source(o.source, {
    format: 'json',
    uri: o.uri,
    fields: { position: 'residueNumber', value: 'confidenceScore' },
    system: o.accession,
    attribution: `AlphaFold DB AF-${o.accession}-F1 v6`,
    licence: 'CC-BY-4.0',
  });
  return b.data.values(o.values, { target: o.target as never, source });
}

/** The prediction's chains: UniProt's entries (their letters and features read by location) and AlphaFold DB's confidence. */
export function prediction(b: SeqSpecBuilder) {
  const chain = (id: 'cdk2' | 'cyclin', accession: 'P24941' | 'P20248', uri: string, confidence: string) => {
    const entry = b.data.source(`${id}-entry`, {
      format: 'json',
      uri,
      ...uniprotTransform({ accession }),
      attribution: 'UniProt Consortium, UniProtKB',
      licence: 'CC-BY-4.0',
    });
    const protein = b.data.sequence(`${id}-protein`, { alphabet: 'protein', source: entry, name: accession });
    const features = b.data.ranges(`${id}-features`, { target: protein, source: entry });
    const plddt = afdbConfidence(b, {
      source: `${id}-plddt`,
      values: `${id}-confidence`,
      accession,
      uri: confidence,
      target: protein,
    });
    /** The data of a layer of one feature type (`Helix`, `Beta strand`). */
    const data = (type: string): string => `${features}${printAddress([{ kind: 'item', params: { category: type } }])}`;
    return { entry, protein, features, plddt, data };
  };
  return {
    cdk2: chain('cdk2', 'P24941', Files.cdk2, Files.cdk2Confidence),
    cyclin: chain('cyclin', 'P20248', Files.cyclin, Files.cyclinConfidence),
  };
}

// ---- the PAE over the complex: synthetic ----

/** How ordered a residue is (1 folded, 0 disordered): cyclin A2's first 172 residues (not in 1FIN's construct) and the chains' ends disordered. */
function order(chain: 'cdk2' | 'cyclin', i: number): number {
  if (chain === 'cdk2') return i < 2 || i >= Lengths.cdk2 - 3 ? 0.3 : 1;
  const first = constructOffset();
  return i < first - 12 ? 0 : i < first ? (i - (first - 12)) / 12 : 1;
}

/**
 * A synthetic predicted aligned error over the complex, AlphaFold-shaped (`[{ predicted_aligned_error: [[…]],
 * max_predicted_aligned_error }]`, Å, whole numbers as AlphaFold DB's v6 files): its rows and columns CDK2's 298
 * residues then cyclin A2's 432. Low within a folded chain, high where either residue is disordered, and across the
 * chains low at 1FIN's interface (its contacts, shifted onto P20248), moderate elsewhere; noise from a seed. Said to be
 * synthetic wherever it is shown.
 */
export function complexPaeDocument(seed = 1995, d: Derived = Data) {
  const n = Lengths.cdk2 + Lengths.cyclin;
  const chainOf = (i: number) => (i < Lengths.cdk2 ? 'cdk2' : 'cyclin');
  const local = (i: number) => (i < Lengths.cdk2 ? i : i - Lengths.cdk2);
  // nearness to the interface: per pair of residues across the chains, the best falloff over 1FIN's contacts
  const shift = constructOffset(d);
  const near = new Float32Array(Lengths.cdk2 * Lengths.cyclin);
  const reach = 10;
  for (const [, a, c] of d.contacts) {
    const b = c + shift;
    for (let i = Math.max(0, a - reach); i <= Math.min(Lengths.cdk2 - 1, a + reach); i++)
      for (let j = Math.max(0, b - reach); j <= Math.min(Lengths.cyclin - 1, b + reach); j++) {
        const f = (1 - Math.abs(i - a) / (reach + 1)) * (1 - Math.abs(j - b) / (reach + 1));
        const k = i * Lengths.cyclin + j;
        if (f > near[k]) near[k] = f;
      }
  }
  const r = new PCG(seed);
  const rows: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < n; j++) {
      const ci = chainOf(i);
      const cj = chainOf(j);
      const o = order(ci, local(i)) * order(cj, local(j));
      const noise = r.float() * 2;
      let v: number;
      if (ci === cj) v = i === j ? 0 : 1.5 + noise + 2.5 * Math.min(1, Math.abs(i - j) / 80) + 26 * (1 - o);
      else {
        const k = ci === 'cdk2' ? local(i) * Lengths.cyclin + local(j) : local(j) * Lengths.cyclin + local(i);
        v = 3 + noise + 8 * (1 - near[k]) + 19 * (1 - o);
      }
      row.push(Math.min(31, Math.round(v)));
    }
    rows.push(row);
  }
  return [{ predicted_aligned_error: rows, max_predicted_aligned_error: 31.75 }];
}

/** The path to the matrix in an AlphaFold PAE document (its first record's member). */
export const PaePath = '0.predicted_aligned_error';

/**
 * The complex's PAE read through the matrix read (the records runtime's, `shape: matrix`), checked against the blocks'
 * total, and cut into the four blocks of chain pairs a grid node each holds (x-major): `[x chain][y chain]`.
 */
export function complexPae(doc: unknown = complexPaeDocument()) {
  const n = Lengths.cdk2 + Lengths.cyclin;
  const r = readMatrix({ value: doc as never }, { records: PaePath } as never, {
    expected: { rows: n, cols: n },
    limits: DefaultOpenLimits,
  });
  if (!r.matrix) throw new Error(`The PAE is not served: ${r.issues.map((i) => i.message).join('; ')}`);
  const m = r.matrix;
  const span = { cdk2: [0, Lengths.cdk2], cyclin: [Lengths.cdk2, n] } as const;
  const cut = (x: 'cdk2' | 'cyclin', y: 'cdk2' | 'cyclin') => {
    const [x0, x1] = span[x];
    const [y0, y1] = span[y];
    // the matrix is row-major by the scored residue (AlphaFold's rows: aligned on i, error at j): x the aligned
    const out = new Float32Array((x1 - x0) * (y1 - y0));
    for (let i = x0; i < x1; i++)
      for (let j = y0; j < y1; j++) out[(i - x0) * (y1 - y0) + (j - y0)] = m.values[i * m.cols + j];
    return out;
  };
  return { matrix: m, block: cut };
}
