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
import type { SchemaItem } from '@molstar/model/model/structure/structure/element/schema';
import type { RowOperator } from '@molstar/model/model/structure/structure/element/schema-selection/selections';

/** An operator for tests which apply rows to a model directly: its instance is '1', as the `instance_id` rows of the fixture expect. */
export const TestOperator: RowOperator = { instanceId: '1', name: '1_555' };
/** The identity operator of a structure made of a model (what its units have) */
export const IdentityOperator: RowOperator = { instanceId: '1_555', name: '1_555' };

export interface Dims {
  nChains: number;
  nResidues: number;
  /** Give every 7th residue the insertion code `A` (off by default to keep the benchmark data fixed). */
  insertionCodes?: boolean;
  /** Make insertion series: the 7th residue of every chain is `<previous auth_seq_id>A` and the 8th is `<previous auth_seq_id>B` (and so on every 7 residues),
   * so that residues share `auth_seq_id` and differ in the insertion code. Replaces `insertionCodes`. */
  insertionSeries?: boolean;
  /** Number the residues of odd chains in descending `label_seq_id` order (so that ranges of `label_seq_id` are not ascending in residue index). */
  descendingSeqIds?: boolean;
  /** Give the last atom of every 5th residue of every 4th chain a different comp_id (`ALT`), so that residues do not have a single comp_id (e.g. microheterogeneity). */
  mixedCompIds?: boolean;
  /** Make `auth_comp_id` differ from `label_comp_id` (prefixed with `X`) in every 3rd chain. */
  differentAuthComp?: boolean;
  /** Let `atom_site.id` wrap around after this many atoms (so that atoms share ids, as in formats with 5 digit serial numbers). */
  atomIdWrap?: number;
  /** Add this many ligands: chains `L0`, `L1`, ... of one HETATM residue (comp_id `LIG`, 4 atoms, auth_seq_id 500 + index) which has no `label_seq_id`, as in real files. */
  ligands?: number;
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
  insertionSeries,
  descendingSeqIds,
  mixedCompIds,
  differentAuthComp,
  atomIdWrap,
  ligands,
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
      // (`A` after a residue, `B` after that: the numbers of an insertion series are the one of the residue before it)
      const authSeq = insertionSeries ? r + 10 - (r % 7 === 0 ? 1 : r % 7 === 1 && r > 1 ? 2 : 0) : r + 10;
      const insCode = insertionSeries
        ? r % 7 === 0
          ? 'A'
          : r % 7 === 1 && r > 1
            ? 'B'
            : '?'
        : insertionCodes && r % 7 === 0
          ? 'A'
          : '?';
      for (let a = 0; a < AtomNames.length; a++) {
        const atomComp = mixedCompIds && c % 4 === 3 && r % 5 === 0 && a === AtomNames.length - 1 ? 'ALT' : comp;
        const x = (c * 10 + a).toFixed(3);
        const y = (r * 3.8).toFixed(3);
        lines.push(
          `ATOM ${atomIdWrap ? ((id++ - 1) % atomIdWrap) + 1 : id++} ${TypeSymbols[a]} ${AtomNames[a]} . ${atomComp} ${asym} ${entity} ${descendingSeqIds && c % 2 === 1 ? NResidues + 1 - r : r} ${insCode} ${x} ${y} 0.000 1 10 ${authSeq} ${differentAuthComp && c % 3 === 1 ? `X${atomComp}` : atomComp} ${auth} ${AtomNames[a]} 1`,
        );
      }
    }
  }
  for (let l = 0; l < (ligands ?? 0); l++) {
    for (let a = 0; a < 4; a++) {
      const serial = atomIdWrap ? ((id++ - 1) % atomIdWrap) + 1 : id++;
      lines.push(
        `HETATM ${serial} ${a % 2 ? 'O' : 'C'} X${a} . LIG L${l} ${NEntities + 1} . ? ${l * 3 + a} -5 0.000 1 10 ${500 + l} LIG L${l} X${a} 1`,
      );
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

export type RowCategory<R = MVSAnnotationRow> = {
  name: string;
  count: number;
  make: (rnd: () => number, i: number) => R;
};

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

export function makeRows<R = MVSAnnotationRow>(
  categories: RowCategory<R>[],
  seed: number,
): { name: string; rows: R[] }[] {
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
  /** A comp_id as it can be in `auth_comp_id` (see `differentAuthComp`) */
  const authComp = (r: () => number) => (r() < 0.5 ? comp(r) : `X${comp(r)}`);
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
    cat('residue_index+label_seq (no chain)', (r) => {
      const i = randInt(r, nChains * nResidues);
      return { residue_index: i, label_seq_id: r() < 0.6 ? 1 + (i % nResidues) : seq(r) };
    }),
    cat('residue_index only, missing', (r) => ({ residue_index: nChains * nResidues + randInt(r, 50) })),
    cat('ins_code+label_comp (no chain, no seq)', (r) => ({ pdbx_PDB_ins_code: 'A', label_comp_id: comp(r) })),
    cat('no insertion code (empty ins code)', (r) => {
      const s = seq(r);
      const pick = r();
      if (pick < 0.25) return { pdbx_PDB_ins_code: '' };
      if (pick < 0.5) return { auth_asym_id: authChain(r), pdbx_PDB_ins_code: '' };
      if (pick < 0.75) return { auth_asym_id: authChain(r), auth_seq_id: s + 10, pdbx_PDB_ins_code: '' };
      return { label_asym_id: chain(r), label_comp_id: comp(r), pdbx_PDB_ins_code: '' };
    }),
    cat('label_comp+auth_comp (no chain)', (r) => ({ label_comp_id: comp(r), auth_comp_id: authComp(r) })),
    cat('auth_comp only (no chain)', (r) => ({ auth_comp_id: authComp(r) })),
    cat('comp+seq (no chain)', (r) => ({ label_comp_id: comp(r), label_seq_id: seq(r) })),
    cat('atom filter only (no chain, no residue)', (r) => ({
      label_atom_id: r() < 0.5 ? 'CA' : 'CG',
      type_symbol: 'C',
    })),
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
    cat('atom_id (repeated ids select all their atoms)', (r) => {
      const i = randInt(r, Math.min(nAtoms, 150));
      const pick = r();
      if (pick < 0.2) return { atom_id: i + 1 };
      if (pick < 0.4) return { atom_id: i + 1, label_asym_id: chain(r) };
      if (pick < 0.6) return { atom_id: i + 1, label_seq_id: seq(r) };
      if (pick < 0.8) return { atom_id: i + 1, atom_index: i };
      return { atom_id: i + 1, atom_index: randInt(r, nAtoms) };
    }),
    cat('instance_id', (r) => ({ instance_id: r() < 0.5 ? '1' : '2', label_asym_id: chain(r), label_seq_id: seq(r) })),
    cat('operator_name', (r) => ({
      operator_name: r() < 0.5 ? '1_555' : '2_555',
      label_asym_id: chain(r),
      label_seq_id: seq(r),
    })),
    cat('instance_id and operator_name', (r) => ({
      instance_id: r() < 0.5 ? '1' : '2',
      operator_name: r() < 0.7 ? '1_555' : '2_555',
      label_asym_id: chain(r),
    })),
  ];
}

/** Rows with ranges of `auth_seq_id` which have insertion codes (and the rows which look like them but are not). Needs `insertionSeries`. */
export function insCodeRangeCategories(
  { nChains, nResidues }: Dims,
  nAtoms: number,
  count = 60,
): RowCategory<SchemaItem>[] {
  const chain = (r: () => number) => chainName(randInt(r, nChains));
  const authChain = (r: () => number) => authChainName(randInt(r, nChains));
  const seq = (r: () => number) => 1 + randInt(r, nResidues);
  /** The `auth_seq_id` of residue `s`, or of a residue in the series after it (their numbers are shared) */
  const auth = (s: number) => s + 10;
  const code = (r: () => number) => ['', 'A', 'B', 'C'][randInt(r, 4)];
  const cat = (name: string, make: RowCategory<SchemaItem>['make']): RowCategory<SchemaItem> => ({ name, count, make });
  const chainCriteria = (r: () => number): SchemaItem => {
    const pick = r();
    return pick < 0.45 ? { label_asym_id: chain(r) } : pick < 0.8 ? { auth_asym_id: authChain(r) } : {};
  };
  /** Bounds near the insertion series (every 7th residue), where codes make a difference */
  const nearSeries = (r: () => number) => {
    const k = 7 * (1 + randInt(r, Math.max(1, Math.floor(nResidues / 7)))) - randInt(r, 2);
    return Math.min(nResidues, Math.max(1, k));
  };
  const range = (r: () => number, startOf: (r: () => number) => number): SchemaItem => {
    const row: SchemaItem = {};
    const a = startOf(r);
    const b = Math.min(nResidues + 3, a + randInt(r, 16));
    const pick = r();
    if (pick < 0.8 || pick >= 0.9) {
      row.beg_auth_seq_id = auth(a);
      if (r() < 0.75) row.beg_auth_ins_code = code(r);
    }
    if (pick < 0.7 || pick >= 0.9) {
      row.end_auth_seq_id = auth(b);
      if (r() < 0.75) row.end_auth_ins_code = code(r);
    }
    return row;
  };
  return [
    cat('chain + range near the series', (r) => ({ ...chainCriteria(r), ...range(r, nearSeries) })),
    cat('chain + range anywhere', (r) => ({ ...chainCriteria(r), ...range(r, seq) })),
    cat('range + residue criteria', (r) => {
      const row = { ...chainCriteria(r), ...range(r, nearSeries) };
      const pick = r();
      if (pick < 0.2) row.label_seq_id = seq(r);
      else if (pick < 0.4) row.pdbx_PDB_ins_code = code(r);
      else if (pick < 0.6) row.residue_index = randInt(r, nChains * nResidues);
      else if (pick < 0.8) row.label_comp_id = CompIds[randInt(r, CompIds.length)];
      else row.beg_label_seq_id = seq(r);
      return row;
    }),
    cat('range + atom criteria', (r) => {
      const row = { ...chainCriteria(r), ...range(r, nearSeries) };
      if (r() < 0.5) row.label_atom_id = 'CA';
      else row.type_symbol = 'c';
      return row;
    }),
    cat('range + atom_id / atom_index (of an atom in the range)', (r) => {
      const c = randInt(r, nChains);
      const res = nearSeries(r);
      const iAtom = (c * nResidues + res - 1) * 6 + randInt(r, 6);
      const row: SchemaItem = {
        label_asym_id: chainName(c),
        beg_auth_seq_id: auth(Math.max(1, res - randInt(r, 3))),
        end_auth_seq_id: auth(res + randInt(r, 3)),
      };
      if (r() < 0.6) row.beg_auth_ins_code = code(r);
      if (r() < 0.6) row.end_auth_ins_code = code(r);
      if (r() < 0.5) row.atom_index = iAtom;
      else row.atom_id = iAtom + 1;
      return row;
    }),
    cat('range + instance', (r) => ({
      ...chainCriteria(r),
      ...range(r, nearSeries),
      instance_id: r() < 0.5 ? '1' : '2',
    })),
    cat('insertion code without its number', (r) => {
      const row: SchemaItem = { ...chainCriteria(r) };
      const pick = r();
      if (pick < 0.4) row.beg_auth_ins_code = code(r);
      else if (pick < 0.8) row.end_auth_ins_code = code(r);
      else {
        row.beg_auth_seq_id = auth(seq(r));
        row.end_auth_ins_code = code(r);
      }
      return row;
    }),
  ];
}
