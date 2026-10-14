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
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { IndicesAndSortings } from '@molstar/model/model/structure/structure/element/schema-selection/indexing';
import { Coordinates, type Frame, Model, StructureElement, Time } from '@molstar/model/model/structure';
import type { MVSAnnotationRow } from '../schemas.js';
import { getAtomRangesForRows } from '@molstar/model/model/structure/structure/element/schema-selection/selections';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import { Structure } from '@molstar/model/model/structure';
import '@molstar/model/script/runtime/query/table'; // (registers the symbols of the queries)
import {
  NEntities,
  TestOperator,
  createModel,
  makeRows,
  mulberry32,
  randInt,
  rowCategories,
} from './annotation-rows-fixture.js';

const NChains = Number(process.env.MVS_BENCH_CHAINS ?? 300);
const NResidues = Number(process.env.MVS_BENCH_RESIDUES ?? 400);
const NRows = Number(process.env.MVS_BENCH_ROWS ?? 1000);
const NRuns = 3;

function allocated() {
  return getHeapStatistics().total_allocated_bytes ?? 0; // (not in the statistics of older versions of Node: allocations then show as 0)
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
  return (annot as any).getRowForEachAtom(model, TestOperator) as {
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

describe.skipIf(process.env.MVS_BENCH !== '1')('annotation rows benchmark', () => {
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

    // rows which select a small part: nothing covers every atom, so every row has to be processed
    // (a set covering everything lets `getRowForEachAtom` stop early)
    const broad = new Set(['comp_id only (all chains)', 'entity only', 'seq only (all chains)']);
    const selective = groups.filter((g) => !broad.has(g.name)).flatMap((g) => g.rows);
    const rndSel = mulberry32(778);
    for (let i = selective.length - 1; i > 0; i--) {
      const j = randInt(rndSel, i + 1);
      [selective[i], selective[j]] = [selective[j], selective[i]];
    }
    run('mixed, selective rows only', selective);

    console.log('\n--- getRowForEachAtom, overlapping rows (entity only; each covers 1/5 of atoms)');
    run(
      'overlap: entity only',
      Array.from({ length: NRows }, (_, i) => ({ label_entity_id: `${1 + (i % NEntities)}` })),
    );

    console.log('\n--- getAtomRangesForRows (union, as used by label-text ElementRangesCache)');
    const model1 = model;
    const rowsArr = mixed;
    const m = measure(() => getAtomRangesForRows(rowsArr, model1, TestOperator, IndicesAndSortings.get(model1)));
    report(
      'union: mixed',
      rowsArr.length,
      m,
      `ranges=${m.result.from.length} checksum=${checksum([m.result.from, m.result.to])}`,
    );

    console.log(
      '\n--- rows -> substructure: `substructureFromRows` vs the expression element-schema (`Schema.toLoci`: expression, compile, query, loci)',
    );
    const structure = Structure.ofModel(model);
    for (const n of [1, 10, 100, 1000, 10000]) {
      const items = (
        n <= mixed.length ? mixed.slice(0, n) : Array.from({ length: n }, (_, i) => mixed[i % mixed.length])
      ) as StructureElement.SchemaItem[];
      const direct = measure(() => substructureFromRows(structure, items));
      const query = measure(() => StructureElement.Schema.toLoci(structure, { items }));
      const size = StructureElement.Loci.size(query.result);
      report(`substructureFromRows (${n} rows)`, n, direct, `elements=${direct.result.elementCount}`);
      report(
        `toLoci, expression     (${n} rows)`,
        n,
        query,
        `elements=${size}  overhead x${(query.ms / Math.max(direct.ms, 0.01)).toFixed(1)}`,
      );
    }

    console.log('\n--- trajectory: 10 frames, 100 rows -> substructure for each frame');
    const nFrames = 10;
    const frames: Frame[] = Array.from({ length: nFrames }, (_, i) => ({
      elementCount: nAtoms,
      time: Time(i, 'step'),
      x: new Float32Array(nAtoms).fill(i),
      y: new Float32Array(nAtoms).fill(i),
      z: new Float32Array(nAtoms).fill(i),
      xyzOrdering: { isIdentity: true },
    }));
    const trajectory = Model.trajectoryFromModelAndCoordinates(
      model,
      Coordinates.create(frames, Time(1, 'step'), Time(0, 'step')),
    );
    const frameModels = Array.from({ length: nFrames }, (_, i) => trajectory.getFrameAtIndex(i) as Model);
    const rows100 = mixed.slice(0, 100);
    IndicesAndSortings.get(model); // (already created above, shared by all frames)
    // per-frame indices (what every frame cost before the indices were shared) + the substructure
    const rebuilt = measure(() => {
      let count = 0;
      for (const f of frameModels) {
        IndicesAndSortings.create(f);
        count += substructureFromRows(Structure.ofModel(f), rows100).elementCount;
      }
      return count;
    });
    const shared = measure(() => {
      let count = 0;
      for (const f of frameModels) count += substructureFromRows(Structure.ofModel(f), rows100).elementCount;
      return count;
    });
    report('frames: indices rebuilt per frame', nFrames, rebuilt, `elements=${rebuilt.result}`);
    report(
      'frames: indices shared',
      nFrames,
      shared,
      `elements=${shared.result}  speedup x${(rebuilt.ms / shared.ms).toFixed(1)}`,
    );
    // the row -> atom mapping of an annotation (colors, tooltips) per frame
    const annotRows = annotationOf(mixed);
    annotRows.getRows();
    const rowMapRebuilt = measure(() => frameModels.map((f) => indexedModelOf(annotationOf(mixed), f).atoms!.length));
    const rowMapShared = measure(() =>
      frameModels.map((f) => (annotRows as any).getIndexedModel(f, TestOperator).atoms.length),
    );
    report('frames: row->atom map rebuilt', nFrames, rowMapRebuilt, '');
    report(
      'frames: row->atom map shared',
      nFrames,
      rowMapShared,
      `speedup x${(rowMapRebuilt.ms / Math.max(rowMapShared.ms, 0.01)).toFixed(0)}`,
    );
  });
});
