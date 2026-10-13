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
import type { Model } from '@molstar/model/model/structure';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { IndicesAndSortings } from '../indexing.js';
import type { MVSAnnotationRow } from '../schemas.js';
import { getAtomRangesForRows, rowsToExpression } from '../selections.js';
import { substructureFromRows } from '../row-structure.js';
import { Structure, StructureSelection } from '@molstar/model/model/structure';
import { StructureQueryHelper } from '@molstar/plugin/state/helpers/structure-query';
import { NEntities, createModel, makeRows, mulberry32, randInt, rowCategories } from './annotation-rows-fixture.js';

const NChains = Number(process.env.MVS_BENCH_CHAINS ?? 300);
const NResidues = Number(process.env.MVS_BENCH_RESIDUES ?? 400);
const NRows = Number(process.env.MVS_BENCH_ROWS ?? 1000);
const NRuns = 3;

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
    const model = await createModel({ nChains: NChains, nResidues: NResidues });
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

    const groups = makeRows(rowCategories({ nChains: NChains, nResidues: NResidues }, nAtoms, NRows), 12345);
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

    console.log(
      '\n--- rows -> substructure: MolScript expression (per-atom queries) vs index-driven (`substructureFromRows`)',
    );
    const structure = Structure.ofModel(model);
    for (const n of [10, 100, 1000]) {
      const rows = mixed.slice(0, n);
      // fewer repetitions for the slow reference
      const ref = measure(() => {
        const { selection } = StructureQueryHelper.createAndRun(structure, rowsToExpression(rows));
        return StructureSelection.unionStructure(selection);
      });
      const fast = measure(() => substructureFromRows(structure, rows));
      report(`expression  (${n} rows)`, n, ref, `elements=${ref.result.elementCount}`);
      report(
        `ranges      (${n} rows)`,
        n,
        fast,
        `elements=${fast.result.elementCount}  speedup x${(ref.ms / fast.ms).toFixed(0)}`,
      );
    }
  });
});
