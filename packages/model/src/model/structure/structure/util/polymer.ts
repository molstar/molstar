/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, type ElementIndex } from '../../../structure.js';
import { Segmentation, SortedArray } from '@molstar/core/data/int';
import { SortedRanges } from '@molstar/core/data/int/sorted-ranges';
import { isNucleic, isProtein } from '../../model/types.js';

export function getAtomicPolymerElements(unit: Unit.Atomic) {
  const indices: ElementIndex[] = [];
  const { elements, model } = unit;
  const { residueAtomSegments } = unit.model.atomicHierarchy;
  const { traceElementIndex } = model.atomicHierarchy.derived.residue;
  const polymerIt = SortedRanges.transientSegments(unit.model.atomicRanges.polymerRanges, elements);
  const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
  while (polymerIt.hasNext) {
    const polymerSegment = polymerIt.move();
    residueIt.setSegment(polymerSegment);
    while (residueIt.hasNext) {
      // the iterator yields only residues with elements in the unit (`start < end`), so every one counts; testing
      // its own element range against the unit's elements was always true (perf-audit-3 05/P25)
      const { index } = residueIt.move();
      const elementIndex = traceElementIndex[index];
      indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex);
    }
  }
  return SortedArray.ofSortedArray<ElementIndex>(indices);
}

export function getCoarsePolymerElements(unit: Unit.Spheres | Unit.Gaussians) {
  const indices: ElementIndex[] = [];
  const { elements, model } = unit;
  const { spheres, gaussians } = model.coarseHierarchy;
  const polymerRanges = Unit.isSpheres(unit) ? spheres.polymerRanges : gaussians.polymerRanges;
  const polymerIt = SortedRanges.transientSegments(polymerRanges, elements);
  while (polymerIt.hasNext) {
    const { start, end } = polymerIt.move();
    for (let i = start; i < end; ++i) {
      indices.push(elements[i]);
    }
  }
  return SortedArray.ofSortedArray<ElementIndex>(indices);
}

//

export function getAtomicGapElements(unit: Unit.Atomic) {
  const indices: ElementIndex[] = [];
  const { elements, model, residueIndex } = unit;
  const { residueAtomSegments } = unit.model.atomicHierarchy;
  const { traceElementIndex } = model.atomicHierarchy.derived.residue;
  const gapIt = SortedRanges.transientSegments(unit.model.atomicRanges.gapRanges, unit.elements);
  while (gapIt.hasNext) {
    const gapSegment = gapIt.move();
    const indexStart = residueIndex[elements[gapSegment.start]];
    const indexEnd = residueIndex[elements[gapSegment.end - 1]];
    const elementIndexStart = traceElementIndex[indexStart];
    const elementIndexEnd = traceElementIndex[indexEnd];
    indices.push(elementIndexStart === -1 ? residueAtomSegments.offsets[indexStart] : elementIndexStart);
    indices.push(elementIndexEnd === -1 ? residueAtomSegments.offsets[indexEnd] : elementIndexEnd);
  }
  return SortedArray.ofSortedArray<ElementIndex>(indices);
}

export function getCoarseGapElements(unit: Unit.Spheres | Unit.Gaussians) {
  const indices: ElementIndex[] = [];
  const { elements, model } = unit;
  const { spheres, gaussians } = model.coarseHierarchy;
  const gapRanges = Unit.isSpheres(unit) ? spheres.gapRanges : gaussians.gapRanges;
  const gapIt = SortedRanges.transientSegments(gapRanges, elements);
  while (gapIt.hasNext) {
    const { start, end } = gapIt.move();
    indices.push(elements[start], elements[end - 1]);
  }
  return SortedArray.ofSortedArray<ElementIndex>(indices);
}

//

export function getNucleotideElements(unit: Unit.Atomic) {
  const indices: ElementIndex[] = [];
  const { elements, model } = unit;
  const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
  const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue;
  const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
  const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
  while (chainIt.hasNext) {
    residueIt.setSegment(chainIt.move());

    while (residueIt.hasNext) {
      const { index } = residueIt.move();

      if (isNucleic(moleculeType[index])) {
        const elementIndex = traceElementIndex[index];
        indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex);
      }
    }
  }
  return SortedArray.ofSortedArray<ElementIndex>(indices);
}

export function getProteinElements(unit: Unit.Atomic) {
  const indices: ElementIndex[] = [];
  const { elements, model } = unit;
  const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
  const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue;
  const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
  const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
  while (chainIt.hasNext) {
    residueIt.setSegment(chainIt.move());

    while (residueIt.hasNext) {
      const { index } = residueIt.move();

      if (isProtein(moleculeType[index])) {
        const elementIndex = traceElementIndex[index];
        indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex);
      }
    }
  }
  return SortedArray.ofSortedArray<ElementIndex>(indices);
}
