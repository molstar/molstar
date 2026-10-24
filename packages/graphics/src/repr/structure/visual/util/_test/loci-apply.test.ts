/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { Interval, OrderedSet, SortedArray } from '@molstar/core/data/int';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { eachInterBond, eachIntraBond, eachStructureGroupsBond } from '../bond.js';
import {
  applyIndexRuns,
  createElementSphereImpostor,
  createStructureElementSphereImpostor,
  eachElement,
} from '../element.js';
import { getGroupLoci, isWholeStructureLoci } from '../loci-index.js';

async function load1cbs() {
  const data = fs.readFileSync(
    path.resolve(import.meta.dirname, '../../../../../../../../data/examples/1cbs_updated.cif'),
    'utf8',
  );
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

/** the intervals `run` applies, as [start, end) pairs */
function collect(run: (apply: (interval: Interval) => boolean) => boolean) {
  const intervals: [number, number][] = [];
  const changed = run((i) => {
    intervals.push([Interval.start(i), Interval.end(i)]);
    return true;
  });
  return { changed, intervals };
}

/** the indices the intervals cover, in order */
function covered(intervals: [number, number][]) {
  const out: number[] = [];
  for (const [s, e] of intervals) for (let i = s; i < e; ++i) out.push(i);
  return out;
}

function element(unit: Unit, indices: number[] | Interval) {
  return {
    unit,
    indices: (Interval.is(indices as Interval)
      ? indices
      : SortedArray.ofSortedArray(indices as number[])) as OrderedSet<StructureElement.UnitIndex>,
  };
}

describe('applyIndexRuns', () => {
  it('applies one interval per run of consecutive indices', () => {
    const { intervals } = collect((apply) =>
      applyIndexRuns(SortedArray.ofSortedArray([1, 2, 3, 7, 9, 10]), 100, apply),
    );
    expect(intervals).toEqual([
      [101, 104],
      [107, 108],
      [109, 111],
    ]);
  });

  it('applies nothing for no indices', () => {
    expect(collect((apply) => applyIndexRuns(SortedArray.ofSortedArray([]), 0, apply))).toEqual({
      changed: false,
      intervals: [],
    });
  });
});

describe('loci application', async () => {
  const model = await load1cbs();
  // two unit cells along a: 8 copies of each unit, one symmetry group per original unit
  const cells = await StructureSymmetry.buildSymmetryRange(model, Vec3.create(0, 0, 0), Vec3.create(1, 0, 0)).run();
  const protein = model.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0) as Unit.Atomic;
  const proteinGroup = model.unitSymmetryGroups[model.unitSymmetryGroupsIndexMap.get(protein.id)];

  it('restricts a loci to the units of a group', () => {
    const group = cells.unitSymmetryGroups.find((g) => g.units.length > 1 && Unit.isAtomic(g.units[0]))!;
    const others = cells.units.filter((u) => !group.unitIndexMap.has(u.id));
    // more elements than the group has units: the indexed path
    const elements = [
      ...others.map((u) => element(u, [0])),
      element(group.units[3], [1]),
      element(group.units[1], [2]),
    ];
    const loci = StructureElement.Loci(cells, elements);
    expect(loci.elements.length).toBeGreaterThan(16);
    const restricted = getGroupLoci(loci, group.units, group.unitIndexMap)!;
    expect(restricted.elements.map((e) => e.unit.id)).toEqual([group.units[3].id, group.units[1].id]);
    expect(restricted.structure).toBe(cells);
    // fewer: the scan
    const small = StructureElement.Loci(cells, [elements[0], elements[elements.length - 1]]);
    expect(getGroupLoci(small, group.units, group.unitIndexMap)!.elements).toEqual([elements[elements.length - 1]]);
    // all of the loci on the group: the loci itself; none: undefined
    const own = StructureElement.Loci(cells, [element(group.units[0], [0])]);
    expect(getGroupLoci(own, group.units, group.unitIndexMap)).toBe(own);
    expect(getGroupLoci(StructureElement.Loci(cells, [elements[0]]), group.units, group.unitIndexMap)).toBeUndefined();
    expect(getGroupLoci(StructureElement.Loci(cells, []), group.units, group.unitIndexMap)).toBeUndefined();
  });

  it('caches whole-structure tests of large loci', () => {
    const whole = StructureElement.Loci.all(cells);
    expect(isWholeStructureLoci(whole)).toBe(true);
    expect(isWholeStructureLoci(whole)).toBe(true);
    const partial = StructureElement.Loci(cells, whole.elements.slice(1));
    expect(isWholeStructureLoci(partial)).toBe(false);
  });

  it('applies the atoms of a residue range as one interval', () => {
    const { intervals } = collect((apply) =>
      eachElement(
        StructureElement.Loci(model, [element(protein, [10, 11, 12, 13, 20, 21])]),
        { structure: model, group: proteinGroup },
        apply,
      ),
    );
    expect(intervals).toEqual([
      [10, 14],
      [20, 22],
    ]);
  });

  describe('intra-unit bonds', () => {
    const { bonds } = protein;
    const groupCount = bonds.edgeCount * 2;
    const unitIdx = proteinGroup.unitIndexMap.get(protein.id);

    /** the directed edges of `indices`, all or (marking) only those with both atoms in `indices`, one by one */
    function expected(indices: number[], isMarking: boolean) {
      const set = new Set(indices);
      const edges: number[] = [];
      for (const v of indices) {
        for (let t = bonds.offset[v]; t < bonds.offset[v + 1]; ++t) {
          if (!isMarking || set.has(bonds.b[t])) edges.push(unitIdx * groupCount + t);
        }
      }
      return edges;
    }

    const cases: [string, number[] | Interval][] = [
      ['whole unit', Interval.ofLength(protein.elements.length)],
      ['atom range', Interval.ofBounds(100, 400)],
      ['scattered atoms', Array.from({ length: 200 }, (_, i) => i * 5)],
      ['one residue', [40, 41, 42, 43, 44, 45, 46, 47]],
      ['nothing', Interval.ofBounds(10, 10)],
    ];

    for (const [name, indices] of cases) {
      for (const isMarking of [true, false]) {
        it(`${name}, ${isMarking ? 'marking' : 'layers'}`, () => {
          const loci = StructureElement.Loci(model, [element(protein, indices)]);
          const list = Interval.is(indices as Interval)
            ? Array.from(
                { length: Interval.size(indices as Interval) },
                (_, i) => Interval.start(indices as Interval) + i,
              )
            : (indices as number[]);
          const want = expected(list, isMarking);

          const intra = collect((apply) =>
            eachIntraBond(loci, { structure: model, group: proteinGroup }, apply, isMarking),
          );
          expect(covered(intra.intervals)).toEqual(want);
          expect(intra.changed).toBe(want.length > 0);
          // runs, not single bonds: no two intervals touch
          for (let i = 1; i < intra.intervals.length; ++i) {
            expect(intra.intervals[i][0]).toBeGreaterThan(intra.intervals[i - 1][1]);
          }

          const o = model.intraUnitBondMapping.unitGroupOffset[model.unitSymmetryGroupsIndexMap.get(protein.id)];
          const groups = collect((apply) => eachStructureGroupsBond(loci, model, apply, isMarking));
          expect(covered(groups.intervals)).toEqual(want.map((t) => t + o));
        });
      }
    }

    it('applies a whole unit as one interval', () => {
      const loci = StructureElement.Loci(model, [element(protein, Interval.ofLength(protein.elements.length))]);
      for (const isMarking of [true, false]) {
        const { intervals } = collect((apply) =>
          eachIntraBond(loci, { structure: model, group: proteinGroup }, apply, isMarking),
        );
        expect(intervals).toEqual([[unitIdx * groupCount, (unitIdx + 1) * groupCount]]);
      }
    });
  });
});

describe('inter-unit bonds of element loci (perf-audit-2 07/P10)', async () => {
  const data = fs.readFileSync(
    path.resolve(import.meta.dirname, '../../../../../../../../data/examples/7qpd.fw2.cif'),
    'utf8',
  );
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const structure = Structure.ofModel(
    (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative,
  );
  const bonds = structure.interUnitBonds;

  /** a loci of every other atom of every unit */
  const loci = StructureElement.Loci(
    structure,
    structure.units.map((unit) => {
      const xs: number[] = [];
      for (let i = 0; i < unit.elements.length; i += 2) xs.push(i);
      return { unit, indices: OrderedSet.ofSortedArray(xs as StructureElement.UnitIndex[]) };
    }),
  );
  const inLoci = (unitId: number, index: number) => index % 2 === 0 && structure.unitMap.has(unitId);

  it('apply the edges a brute-force scan finds, without scanning the connected atoms', () => {
    expect(bonds.edgeCount).toBeGreaterThan(0);
    const spies = structure.units.flatMap((u) =>
      bonds.getConnectedUnits(u.id).map((b) => vi.spyOn(b.connectedIndices as number[], 'includes')),
    );
    for (const isMarking of [false, true]) {
      const expected: number[] = [];
      for (const e of bonds.edges) {
        if (inLoci(e.unitA, e.indexA) && (!isMarking || inLoci(e.unitB, e.indexB))) {
          expected.push(bonds.getEdgeIndex(e.indexA, e.unitA, e.indexB, e.unitB));
        }
      }
      const { intervals } = collect((apply) => eachInterBond(loci, structure, apply, isMarking));
      expect(intervals.map((i) => i[0]).sort((a, b) => a - b)).toEqual(expected.sort((a, b) => a - b));
    }
    for (const s of spies) {
      expect(s).not.toHaveBeenCalled();
      s.mockRestore();
    }
  });
});

describe('sphere impostors (perf-audit-2 07/P18)', async () => {
  const structure = await load1cbs();
  const props = { sizeFactor: 1, ignoreHydrogens: false, ignoreHydrogensVariant: 'all' as const, traceOnly: false };
  let sizes = 0;
  const theme = {
    color: {} as any,
    size: { granularity: 'group', size: () => (sizes++, 1.5) } as any,
  };

  it('evaluate the size theme only for a new bounding sphere', () => {
    const unit = structure.units[0];
    sizes = 0;
    const first = createElementSphereImpostor({} as any, unit, structure, theme, props);
    expect(sizes).toBe(unit.elements.length);
    const radius = first.boundingSphere.radius;
    sizes = 0;
    // a rebuild over the old spheres keeps their bounding sphere
    const second = createElementSphereImpostor({} as any, unit, structure, theme, props, first);
    expect(sizes).toBe(0);
    expect(second.boundingSphere.radius).toBe(radius);

    sizes = 0;
    const all = createStructureElementSphereImpostor({} as any, structure, theme, props);
    expect(sizes).toBe(structure.elementCount);
    sizes = 0;
    createStructureElementSphereImpostor({} as any, structure, theme, props, all);
    expect(sizes).toBe(0);
  });
});
