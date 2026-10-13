/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Synthetic multi-chain structure and a mix of annotation rows, shared by the annotation-rows tests and benchmark. */

import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import type { Model } from '@molstar/model/model/structure';
import type { MVSAnnotationRow } from '../schemas.js';

export interface Dims {
  nChains: number;
  nResidues: number;
  /** Give every 7th residue the insertion code `A` (off by default to keep the benchmark data fixed). */
  insertionCodes?: boolean;
  /** Number the residues of odd chains in descending `label_seq_id` order (so that ranges of `label_seq_id` are not ascending in residue index). */
  descendingSeqIds?: boolean;
  /** Give the last atom of every 5th residue of every 4th chain a different comp_id (`ALT`), so that residues do not have a single comp_id (e.g. microheterogeneity). */
  mixedCompIds?: boolean;
}

export const CompIds = ['ALA', 'GLY', 'LEU', 'SER', 'VAL'];
const AtomNames = ['N', 'CA', 'C', 'O', 'CB', 'CG'];
const TypeSymbols = ['N', 'C', 'C', 'O', 'C', 'C'];
export const NEntities = 5;

/** `0 -> A, 25 -> Z, 26 -> AA, ...` */
export function chainName(i: number): string {
  let s = '';
  i += 1;
  while (i > 0) {
    s = String.fromCharCode(65 + ((i - 1) % 26)) + s;
    i = Math.floor((i - 1) / 26);
  }
  return s;
}
/** Every third chain shares its auth_asym_id with the previous one (as in assemblies / split chains). */
export function authChainName(i: number): string {
  return `x${i - (i % 3 === 2 ? 1 : 0)}`;
}

export function mulberry32(seed: number) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function createMmcifText({
  nChains: NChains,
  nResidues: NResidues,
  insertionCodes,
  descendingSeqIds,
  mixedCompIds,
}: Dims): string {
  const lines: string[] = [
    'data_BENCH',
    'loop_',
    '_atom_site.group_PDB',
    '_atom_site.id',
    '_atom_site.type_symbol',
    '_atom_site.label_atom_id',
    '_atom_site.label_alt_id',
    '_atom_site.label_comp_id',
    '_atom_site.label_asym_id',
    '_atom_site.label_entity_id',
    '_atom_site.label_seq_id',
    '_atom_site.pdbx_PDB_ins_code',
    '_atom_site.Cartn_x',
    '_atom_site.Cartn_y',
    '_atom_site.Cartn_z',
    '_atom_site.occupancy',
    '_atom_site.B_iso_or_equiv',
    '_atom_site.auth_seq_id',
    '_atom_site.auth_comp_id',
    '_atom_site.auth_asym_id',
    '_atom_site.auth_atom_id',
    '_atom_site.pdbx_PDB_model_num',
  ];
  let id = 1;
  for (let c = 0; c < NChains; c++) {
    const asym = chainName(c);
    const auth = authChainName(c);
    const entity = (c % NEntities) + 1;
    for (let r = 1; r <= NResidues; r++) {
      const comp = CompIds[(r + c) % CompIds.length];
      for (let a = 0; a < AtomNames.length; a++) {
        const atomComp = mixedCompIds && c % 4 === 3 && r % 5 === 0 && a === AtomNames.length - 1 ? 'ALT' : comp;
        const x = (c * 10 + a).toFixed(3);
        const y = (r * 3.8).toFixed(3);
        lines.push(
          `ATOM ${id++} ${TypeSymbols[a]} ${AtomNames[a]} . ${atomComp} ${asym} ${entity} ${descendingSeqIds && c % 2 === 1 ? NResidues + 1 - r : r} ${insertionCodes && r % 7 === 0 ? 'A' : '?'} ${x} ${y} 0.000 1 10 ${r + 10} ${atomComp} ${auth} ${AtomNames[a]} 1`,
        );
      }
    }
  }
  return lines.join('\n') + '\n';
}

export async function createModel(dims: Dims): Promise<Model> {
  const parsed = await CIF.parseText(createMmcifText(dims)).run();
  if (parsed.isError) throw new Error(parsed.message);
  const file = parsed.result;
  const trajectory = await trajectoryFromMmCIF(file.blocks[0], file).run();
  return trajectory.representative;
}

export type RowCategory = { name: string; count: number; make: (rnd: () => number, i: number) => MVSAnnotationRow };

export function randInt(rnd: () => number, n: number) {
  return Math.floor(rnd() * n);
}

/** Mix of selector kinds; `count`s are scaled to `nRows` (they are per 1000 rows). */
export function rowCategories(
  { nChains: NChains, nResidues: NResidues }: Dims,
  nAtoms: number,
  NRows: number,
): RowCategory[] {
  const chain = (rnd: () => number) => chainName(randInt(rnd, NChains));
  const seq = (rnd: () => number) => 1 + randInt(rnd, NResidues);
  const span = (rnd: () => number) => {
    const beg = seq(rnd);
    return { beg_label_seq_id: beg, end_label_seq_id: Math.min(NResidues, beg + 5 + randInt(rnd, 35)) };
  };
  return [
    { name: 'residue (asym+seq)', count: 250, make: (r) => ({ label_asym_id: chain(r), label_seq_id: seq(r) }) },
    { name: 'range (asym+beg/end)', count: 150, make: (r) => ({ label_asym_id: chain(r), ...span(r) }) },
    {
      name: 'atom (asym+seq+atom_id)',
      count: 100,
      make: (r) => ({ label_asym_id: chain(r), label_seq_id: seq(r), label_atom_id: 'CA' }),
    },
    {
      name: 'residue_index only',
      count: 100,
      make: (r) => ({ residue_index: randInt(r, NChains * NResidues) }),
    },
    { name: 'comp_id only (all chains)', count: 10, make: (r) => ({ label_comp_id: CompIds[randInt(r, 5)] }) },
    {
      name: 'asym+comp_id',
      count: 40,
      make: (r) => ({ label_asym_id: chain(r), label_comp_id: CompIds[randInt(r, 5)] }),
    },
    {
      name: 'auth (auth_asym+auth_seq)',
      count: 100,
      make: (r) => ({ auth_asym_id: authChainName(randInt(r, NChains)), auth_seq_id: seq(r) + 10 }),
    },
    { name: 'entity only', count: 10, make: (r) => ({ label_entity_id: `${1 + randInt(r, NEntities)}` }) },
    { name: 'chain only', count: 30, make: (r) => ({ label_asym_id: chain(r) }) },
    { name: 'atom_index', count: 100, make: (r) => ({ atom_index: randInt(r, nAtoms) }) },
    {
      name: 'asym+atom_id (whole chain CA)',
      count: 40,
      make: (r) => ({ label_asym_id: chain(r), label_atom_id: 'CA' }),
    },
    { name: 'seq only (all chains)', count: 20, make: (r) => ({ label_seq_id: seq(r) }) },
    {
      name: 'asym+range+type_symbol',
      count: 50,
      make: (r) => ({ label_asym_id: chain(r), ...span(r), type_symbol: 'N' }),
    },
  ].map((c) => ({ ...c, count: Math.max(1, Math.round((c.count * NRows) / 1000)) }));
}

export function makeRows(categories: RowCategory[], seed: number): { name: string; rows: MVSAnnotationRow[] }[] {
  const rnd = mulberry32(seed);
  return categories.map((c) => ({ name: c.name, rows: Array.from({ length: c.count }, (_, i) => c.make(rnd, i)) }));
}

/** Rows combining several criteria (some contradictory on purpose, so that they select nothing), which the plain mix above does not cover:
 * every criterion which is not the one used to find candidates has to be checked on its own. Needs `insertionCodes`. */
export function combinedCategories({ nChains, nResidues }: Dims, nAtoms: number, count = 40): RowCategory[] {
  const chain = (r: () => number) => chainName(randInt(r, nChains));
  const authChain = (r: () => number) => authChainName(randInt(r, nChains));
  const seq = (r: () => number) => 1 + randInt(r, nResidues);
  /** Matches `label_seq_id` of the same residue most of the time, so that combined rows select something. */
  const authOf = (r: () => number, labelSeq: number) => (r() < 0.7 ? labelSeq + 10 : seq(r) + 10);
  const comp = (r: () => number) => (r() < 0.15 ? 'ALT' : CompIds[randInt(r, CompIds.length)]);
  const entity = (r: () => number) => `${1 + randInt(r, NEntities)}`;
  const span = (r: () => number) => {
    const a = seq(r);
    const b = a + randInt(r, 20);
    return r() < 0.5 ? { lo: a, hi: b } : { lo: a, hi: undefined as number | undefined };
  };
  const cat = (name: string, make: RowCategory['make']): RowCategory => ({ name, count, make });
  return [
    cat('label_asym+label_seq+auth_seq', (r) => {
      const s = seq(r);
      return { label_asym_id: chain(r), label_seq_id: s, auth_seq_id: authOf(r, s) };
    }),
    cat('auth_asym+auth_seq+label_seq', (r) => {
      const s = seq(r);
      return { auth_asym_id: authChain(r), auth_seq_id: authOf(r, s), label_seq_id: s };
    }),
    cat('residue_index+label_asym', (r) => {
      const i = randInt(r, nChains * nResidues);
      const c = r() < 0.7 ? chainName(Math.floor(i / nResidues)) : chain(r);
      return { residue_index: i, label_asym_id: c };
    }),
    cat('residue_index+auth_seq+label_comp', (r) => ({
      residue_index: randInt(r, nChains * nResidues),
      auth_seq_id: seq(r) + 10,
      label_comp_id: comp(r),
    })),
    cat('label_seq+label_comp', (r) => ({ label_asym_id: chain(r), label_seq_id: seq(r), label_comp_id: comp(r) })),
    cat('label_range+auth_range', (r) => {
      const l = span(r);
      const a = span(r);
      return {
        label_asym_id: chain(r),
        beg_label_seq_id: l.lo,
        end_label_seq_id: l.hi,
        beg_auth_seq_id: a.lo + 10,
        end_auth_seq_id: a.hi === undefined ? undefined : a.hi + 10,
      };
    }),
    cat('auth_range+label_comp', (r) => {
      const a = span(r);
      return {
        auth_asym_id: authChain(r),
        beg_auth_seq_id: a.lo + 10,
        end_auth_seq_id: a.hi === undefined ? undefined : a.hi + 10,
        label_comp_id: comp(r),
      };
    }),
    cat('auth_seq+label_range', (r) => {
      const l = span(r);
      return { label_asym_id: chain(r), auth_seq_id: seq(r) + 10, beg_label_seq_id: l.lo, end_label_seq_id: l.hi };
    }),
    cat('residue_index+label_range', (r) => {
      const l = span(r);
      return { residue_index: randInt(r, nChains * nResidues), beg_label_seq_id: l.lo, end_label_seq_id: l.hi };
    }),
    cat('ins_code+label_range', (r) => {
      const l = span(r);
      return { pdbx_PDB_ins_code: 'A', beg_label_seq_id: l.lo, end_label_seq_id: l.hi };
    }),
    cat('label_range, one-sided', (r) => {
      const s = seq(r);
      return r() < 0.5
        ? { label_asym_id: chain(r), beg_label_seq_id: s }
        : { label_asym_id: chain(r), end_label_seq_id: s };
    }),
    cat('ins_code+label_asym(+seq)', (r) => {
      const row: MVSAnnotationRow = { pdbx_PDB_ins_code: 'A', label_asym_id: chain(r) };
      if (r() < 0.5) row.label_seq_id = 7 * (1 + randInt(r, Math.max(1, Math.floor(nResidues / 7))));
      return row;
    }),
    cat('ins_code+auth_seq / ins_code only', (r) =>
      r() < 0.5
        ? { pdbx_PDB_ins_code: 'A' }
        : {
            pdbx_PDB_ins_code: r() < 0.8 ? 'A' : 'B',
            auth_seq_id: 7 * (1 + randInt(r, Math.floor(nResidues / 7))) + 10,
          },
    ),
    cat('auth_comp+label_asym / auth_comp only', (r) =>
      r() < 0.5 ? { auth_comp_id: comp(r), label_asym_id: chain(r) } : { auth_comp_id: comp(r) },
    ),
    cat('label_comp+auth_comp', (r) => ({ label_asym_id: chain(r), label_comp_id: comp(r), auth_comp_id: comp(r) })),
    cat('label_comp+label_atom', (r) => ({ label_comp_id: comp(r), label_atom_id: 'CA', label_asym_id: chain(r) })),
    cat('auth_atom+type_symbol+label_asym', (r) => ({
      auth_atom_id: r() < 0.5 ? 'CB' : 'N',
      type_symbol: r() < 0.5 ? 'c' : 'N',
      label_asym_id: chain(r),
    })),
    cat('label_asym+auth_asym', (r) => ({ label_asym_id: chain(r), auth_asym_id: authChain(r) })),
    cat('label_entity+label_asym', (r) => ({ label_entity_id: entity(r), label_asym_id: chain(r) })),
    cat('label_entity+auth_asym+label_asym', (r) => ({
      label_entity_id: entity(r),
      auth_asym_id: authChain(r),
      label_asym_id: chain(r),
    })),
    cat('atom_id+label_asym / atom_index+label_seq', (r) => {
      const i = randInt(r, nAtoms);
      return r() < 0.5 ? { atom_id: i + 1, label_asym_id: chain(r) } : { atom_index: i, label_seq_id: seq(r) };
    }),
    cat('instance_id', (r) => ({ instance_id: r() < 0.5 ? '1' : '2', label_asym_id: chain(r), label_seq_id: seq(r) })),
  ];
}
