/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { Coordinates, type Frame, Model, Structure, StructureElement, Time } from '@molstar/model/model/structure';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { HierarchyScope, IndicesAndSortings } from '../indexing.js';
import { substructureFromRows } from '../row-structure.js';
import { createModel } from './annotation-rows-fixture.js';

/** A coordinate trajectory of `nFrames` frames: copies of `model` with other coordinates, as created for XTC, DCD, ... */
function trajectoryOf(model: Model, nFrames: number): Model[] {
  const n = model.atomicHierarchy.atoms._rowCount;
  const frames: Frame[] = Array.from({ length: nFrames }, (_, i) => ({
    elementCount: n,
    time: Time(i, 'step'),
    x: new Float32Array(n).fill(i),
    y: new Float32Array(n).fill(2 * i),
    z: new Float32Array(n).fill(3 * i),
    xyzOrdering: { isIdentity: true },
  }));
  const trajectory = Model.trajectoryFromModelAndCoordinates(
    model,
    Coordinates.create(frames, Time(1, 'step'), Time(0, 'step')),
  );
  return Array.from({ length: trajectory.frameCount }, (_, i) => trajectory.getFrameAtIndex(i) as Model);
}

describe('indices are shared by the frames of a trajectory', async () => {
  const model = await createModel({ nChains: 5, nResidues: 20 });
  const frames = trajectoryOf(model, 4);

  it('the frames are different models with the same hierarchies', () => {
    expect(new Set(frames.map((f) => f.id)).size).toBe(4);
    expect(frames.every((f) => Model.areHierarchiesEqual(f, frames[0]))).toBe(true);
    expect(frames[1].atomicConformation).not.toBe(frames[0].atomicConformation);
    expect(frames[1].atomicConformation.x).not.toEqual(frames[0].atomicConformation.x);
  });

  it('HierarchyScope and IndicesAndSortings are the same object for all frames (and the original model)', () => {
    const scope = HierarchyScope.get(frames[0]);
    for (const f of frames) expect(HierarchyScope.get(f)).toBe(scope);
    expect(HierarchyScope.get(model)).toBe(scope);
    const indices = IndicesAndSortings.get(frames[0]);
    for (const f of frames) expect(IndicesAndSortings.get(f)).toBe(indices);
    expect(IndicesAndSortings.get(model)).toBe(indices);
  });

  it('indices are created once for all frames', () => {
    const fresh = trajectoryOf(model, 3).map((f) => ({
      ...f,
      id: `${f.id}-copy`,
      _dynamicPropertyData: Object.create(null),
    }));
    const create = vi.spyOn(IndicesAndSortings, 'create');
    try {
      for (const f of fresh) IndicesAndSortings.get(f as Model);
      // the original model created them (shared hierarchy), nothing to create for new frames
      expect(create).not.toHaveBeenCalled();
    } finally {
      create.mockRestore();
    }
  });

  it('models which only look alike do not share: other hierarchy, or the same hierarchy but other atom ids', async () => {
    const other = await createModel({ nChains: 5, nResidues: 20 }); // same content, own hierarchy objects
    expect(HierarchyScope.get(other)).not.toBe(HierarchyScope.get(model));
    expect(IndicesAndSortings.get(other)).not.toBe(IndicesAndSortings.get(model));

    const otherAtomIds = Column.ofIntArray(
      Array.from({ length: model.atomicHierarchy.atoms._rowCount }, (_, i) => i + 1000),
    );
    const renumbered = { ...model, atomicConformation: { ...model.atomicConformation, atomId: otherAtomIds } } as Model;
    expect(Model.areHierarchiesEqual(renumbered, model)).toBe(true);
    expect(HierarchyScope.get(renumbered)).not.toBe(HierarchyScope.get(model));
    // both stay valid, looking one up does not evict the other
    expect(HierarchyScope.get(model)).toBe(HierarchyScope.get(frames[0]));
    expect(IndicesAndSortings.get(renumbered).atomic!.atomsById.get(1000)).toBe(0);
    expect(IndicesAndSortings.get(model).atomic!.atomsById.get(1000)).toBeUndefined();
  });
});

describe('annotation rows are resolved once for all frames of a trajectory', async () => {
  const model = await createModel({ nChains: 5, nResidues: 20 });
  const frames = trajectoryOf(model, 3);
  const rows = [
    { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 5 },
    { label_asym_id: 'B', label_seq_id: 7, label_atom_id: 'CA' },
  ];
  const valueAt = (annot: MVSAnnotation, frame: Model, iAtom: number, field: string) => {
    const structure = Structure.ofModel(frame);
    const unit = structure.units[0];
    const loc = StructureElement.Location.create(structure, unit, iAtom as any);
    return annot.getValueForLocation(loc, field);
  };

  it('the row of each atom is computed once and is the same in every frame', () => {
    const annot = new MVSAnnotation(
      { format: 'json', data: rows.map((r, i) => ({ ...r, tag: `row${i}` })) as any },
      'all_atomic',
      {},
    );
    const compute = vi.spyOn(annot as any, 'getRowForEachAtom');
    // chain A (atoms 0..119): residues 1..5 are atoms 0..29; chain B's residue 7, atom CA is atom 120 + 6*6 + 1
    const probes = [
      [0, 'row0'],
      [29, 'row0'],
      [30, undefined],
      [120 + 6 * 6 + 1, 'row1'],
    ] as const;
    for (const f of frames) {
      for (const [iAtom, expected] of probes) expect(valueAt(annot, f, iAtom, 'tag') ?? undefined).toBe(expected);
    }
    expect(compute).toHaveBeenCalledTimes(1);
  });

  it('with instance ids, the mapping is still computed once per instance and frame-independent', () => {
    const annot = new MVSAnnotation(
      { format: 'json', data: rows.map((r, i) => ({ ...r, instance_id: '1_555', tag: `row${i}` })) as any },
      'all_atomic',
      {},
    );
    const compute = vi.spyOn(annot as any, 'getRowForEachAtom');
    for (const f of frames) expect(valueAt(annot, f, 0, 'tag')).toBe('row0');
    expect(compute).toHaveBeenCalledTimes(1);
  });

  it('substructureFromRows gives the same selection in every frame, with the coordinates of that frame', () => {
    const expected = Array.from(substructureFromRows(Structure.ofModel(frames[0]), rows).units[0].elements);
    expect(expected.length).toBeGreaterThan(0);
    frames.forEach((f, i) => {
      const sub = substructureFromRows(Structure.ofModel(f), rows);
      expect(Array.from(sub.units[0].elements)).toEqual(expected);
      const unit = sub.units[0];
      expect(unit.conformation.x(unit.elements[0])).toBeCloseTo(i, 5); // frame i has x = i
    });
  });
});
