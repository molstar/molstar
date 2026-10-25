/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Interval, OrderedSet, Segmentation, SortedArray } from '@molstar/core/data/int';
import { SortedRanges } from '@molstar/core/data/int/sorted-ranges';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import type { ElementIndex } from '../../../model/indexing.js';
import { Structure } from '../../structure.js';
import { Unit } from '../../unit.js';
import { getAtomicPolymerElements } from '../polymer.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

/** The former version, with the per-residue intersection test */
function reference(unit: Unit.Atomic) {
  const indices: ElementIndex[] = [];
  const { elements, model } = unit;
  const { residueAtomSegments } = model.atomicHierarchy;
  const { traceElementIndex } = model.atomicHierarchy.derived.residue;
  const polymerIt = SortedRanges.transientSegments(model.atomicRanges.polymerRanges, elements);
  const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
  while (polymerIt.hasNext) {
    residueIt.setSegment(polymerIt.move());
    while (residueIt.hasNext) {
      const { start, end, index } = residueIt.move();
      if (OrderedSet.areIntersecting(Interval.ofRange(elements[start], elements[end - 1]), elements)) {
        const elementIndex = traceElementIndex[index];
        indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex);
      }
    }
  }
  return SortedArray.ofSortedArray<ElementIndex>(indices);
}

describe('getAtomicPolymerElements', () => {
  it('gives the trace elements of every residue with atoms in the unit, as with the dropped test (perf-audit-3 05/P25)', async () => {
    const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const root = Structure.ofModel(model);
    const unit = root.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0) as Unit.Atomic;
    const all = Array.from(unit.elements);

    // whole, every k-th atom (residues with holes, residues left out), one atom, a run inside one residue
    const subsets: ElementIndex[][] = [all, ...[2, 3, 7, 11, 40].map((k) => all.filter((_, i) => i % k === 0))];
    subsets.push([all[100]], all.slice(3, 6));
    let traces = 0;
    for (const s of subsets) {
      const child = unit.getChild(SortedArray.ofSortedArray(s)) as Unit.Atomic;
      const got = getAtomicPolymerElements(child);
      expect(Array.from(got)).toEqual(Array.from(reference(child)));
      traces += got.length;
    }
    expect(traces).toBeGreaterThan(unit.polymerElements.length);
  });
});
