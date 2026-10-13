/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Benchmark for `MVSAnnotation.getRowForEachAtom` (rows -> per-atom row index) and `getAtomRangesForRows` (union of
 * rows) on a synthetic multi-chain structure. Skipped unless `MVS_BENCH=1`:
 *
 *   MVS_BENCH=1 pnpm test:unit packages/mvs/runtime/src/helpers/_test/annotation-rows.bench.test.ts
 *
 * Size knobs: `MVS_BENCH_CHAINS` (default 300), `MVS_BENCH_RESIDUES` (default 400), `MVS_BENCH_ROWS` (default 1000).
 * Each case prints best-of-N time, bytes allocated by one run (`v8` `total_allocated_bytes`), the number of atoms
 * that received a row, and a checksum of the per-atom result (must stay the same across refactors).
 */

import { getHeapStatistics } from 'node:v8';
import { describe, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import type { Model } from '@molstar/model/model/structure';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { IndicesAndSortings } from '../indexing.js';
import type { MVSAnnotationRow } from '../schemas.js';
import { getAtomRangesForRows } from '../selections.js';

const NChains = Number(process.env.MVS_BENCH_CHAINS ?? 300);
const NResidues = Number(process.env.MVS_BENCH_RESIDUES ?? 400);
const NRows = Number(process.env.MVS_BENCH_ROWS ?? 1000);
const NRuns = 3;

const CompIds = ['ALA', 'GLY', 'LEU', 'SER', 'VAL'];
const AtomNames = ['N', 'CA', 'C', 'O', 'CB', 'CG'];
const TypeSymbols = ['N', 'C', 'C', 'O', 'C', 'C'];
const NEntities = 5;

/** `0 -> A, 25 -> Z, 26 -> AA, ...` */
function chainName(i: number): string {
  let s = '';
  i += 1;
  while (i > 0) {
    s = String.fromCharCode(65 + ((i - 1) % 26)) + s;
    i = Math.floor((i - 1) / 26);
  }
  return s;
}
/** Every third chain shares its auth_asym_id with the previous one (as in assemblies / split chains). */
function authChainName(i: number): string {
  return `x${i - (i % 3 === 2 ? 1 : 0)}`;
}

function mulberry32(seed: number) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function createMmcifText(): string {
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
        const x = (c * 10 + a).toFixed(3);
        const y = (r * 3.8).toFixed(3);
        lines.push(
          `ATOM ${id++} ${TypeSymbols[a]} ${AtomNames[a]} . ${comp} ${asym} ${entity} ${r} ? ${x} ${y} 0.000 1 10 ${r + 10} ${comp} ${auth} ${AtomNames[a]} 1`,
        );
      }
    }
  }
  return lines.join('\n') + '\n';
}

async function createModel(): Promise<Model> {
  const parsed = await CIF.parseText(createMmcifText()).run();
  if (parsed.isError) throw new Error(parsed.message);
  const file = parsed.result;
  const trajectory = await trajectoryFromMmCIF(file.blocks[0], file).run();
  return trajectory.representative;
}

type RowCategory = { name: string; count: number; make: (rnd: () => number, i: number) => MVSAnnotationRow };

function randInt(rnd: () => number, n: number) {
  return Math.floor(rnd() * n);
}

/** Mix of selector kinds; `count`s are scaled to `NRows`. */
function rowCategories(nAtoms: number): RowCategory[] {
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

function makeRows(categories: RowCategory[], seed: number): { name: string; rows: MVSAnnotationRow[] }[] {
  const rnd = mulberry32(seed);
  return categories.map((c) => ({ name: c.name, rows: Array.from({ length: c.count }, (_, i) => c.make(rnd, i)) }));
}

function allocated() {
  return getHeapStatistics().total_allocated_bytes;
}

function measure<T>(f: () => T): { result: T; ms: number; mb: number } {
  let best = Infinity;
  let mb = 0;
  let result!: T;
  for (let i = 0; i < NRuns; i++) {
    const a0 = allocated();
    const t0 = performance.now();
    result = f();
    const dt = performance.now() - t0;
    if (i === 0) mb = (allocated() - a0) / 1e6;
    best = Math.min(best, dt);
  }
  return { result, ms: best, mb };
}

function checksum(arrays: (ArrayLike<number> | null)[]): string {
  let h = 2166136261;
  for (const a of arrays) {
    if (!a) {
      h = Math.imul(h ^ 0xffff, 16777619);
      continue;
    }
    for (let i = 0; i < a.length; i++) h = Math.imul(h ^ (a[i] + 1), 16777619);
  }
  return (h >>> 0).toString(16);
}

function annotationOf(rows: MVSAnnotationRow[]) {
  return new MVSAnnotation({ format: 'json', data: rows as any }, 'all_atomic', {});
}

function indexedModelOf(annot: MVSAnnotation, model: Model) {
  return (annot as any).getRowForEachAtom(model, '1') as {
    atoms: ArrayLike<number> | null;
    spheres: ArrayLike<number> | null;
    gaussians: ArrayLike<number> | null;
  };
}

function report(label: string, nRows: number, m: { ms: number; mb: number }, extra: string) {
  console.log(
    `${label.padEnd(34)} rows=${String(nRows).padStart(5)}  ${m.ms.toFixed(1).padStart(9)} ms  ${m.mb.toFixed(1).padStart(8)} MB alloc  ${extra}`,
  );
}

describe.skipIf(!process.env.MVS_BENCH)('annotation rows benchmark', () => {
  it('getRowForEachAtom / getAtomRangesForRows', { timeout: 600_000 }, async () => {
    let t0 = performance.now();
    const model = await createModel();
    const nAtoms = model.atomicHierarchy.atoms._rowCount;
    console.log(
      `structure: ${NChains} chains x ${NResidues} residues, ${nAtoms} atoms (created in ${(performance.now() - t0).toFixed(0)} ms)`,
    );
    t0 = performance.now();
    const a0 = allocated();
    IndicesAndSortings.get(model);
    console.log(
      `IndicesAndSortings.get: ${(performance.now() - t0).toFixed(0)} ms, ${((allocated() - a0) / 1e6).toFixed(1)} MB alloc (one-off, excluded below)\n`,
    );

    const groups = makeRows(rowCategories(nAtoms), 12345);
    const run = (label: string, rows: MVSAnnotationRow[]) => {
      const annot = annotationOf(rows);
      annot.getRows(); // parse rows outside of the measured part
      const m = measure(() => indexedModelOf(annot, model));
      const atoms = m.result.atoms;
      let assigned = 0;
      if (atoms) for (let i = 0; i < atoms.length; i++) if (atoms[i] >= 0) assigned++;
      report(label, rows.length, m, `assigned=${assigned} checksum=${checksum([atoms])}`);
    };

    console.log('--- getRowForEachAtom, per category');
    for (const g of groups) run(g.name, g.rows);

    console.log('\n--- getRowForEachAtom, mixed');
    const mixed = groups.flatMap((g) => g.rows);
    // interleave deterministically so that categories alternate (as in real annotation files)
    const rnd = mulberry32(777);
    for (let i = mixed.length - 1; i > 0; i--) {
      const j = randInt(rnd, i + 1);
      [mixed[i], mixed[j]] = [mixed[j], mixed[i]];
    }
    run('mixed', mixed);

    console.log('\n--- getRowForEachAtom, overlapping rows (entity only; each covers 1/5 of atoms)');
    run(
      'overlap: entity only',
      Array.from({ length: NRows }, (_, i) => ({ label_entity_id: `${1 + (i % NEntities)}` })),
    );

    console.log('\n--- getAtomRangesForRows (union, as used by label-text ElementRangesCache)');
    const model1 = model;
    const rowsArr = mixed;
    const m = measure(() => getAtomRangesForRows(rowsArr, model1, '1', IndicesAndSortings.get(model1)));
    report(
      'union: mixed',
      rowsArr.length,
      m,
      `ranges=${m.result.from.length} checksum=${checksum([m.result.from, m.result.to])}`,
    );
  });
});
