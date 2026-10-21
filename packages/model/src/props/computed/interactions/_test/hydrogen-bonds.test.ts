/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SortedArray } from '@molstar/core/data/int';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { Structure, type StructureElement, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { eachIntraBondedAtom, typeSymbol } from '../../chemistry/util.js';
import { FeatureType, InteractionType } from '../common.js';
import { computeInteractions } from '../interactions.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');
const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };

/** Donor–acceptor element pairs of the intra- and inter-unit H-bonds */
async function hydrogenBonds(structure: Structure) {
  const { contacts, unitsContacts, unitsFeatures } = await computeInteractions(ctx, structure, {});
  const ret = new Set<string>();
  const element = (unit: Unit, feature: number) => {
    const { offsets, members } = unitsFeatures.get(unit.id);
    return unit.elements[members[offsets[feature]]];
  };
  for (const e of contacts.edges) {
    if (e.props.type !== InteractionType.HydrogenBond) continue;
    ret.add(
      `${element(structure.unitMap.get(e.unitA), e.indexA)}-${element(structure.unitMap.get(e.unitB), e.indexB)}`,
    );
  }
  for (const unit of structure.units) {
    const c = unitsContacts.get(unit.id);
    for (let i = 0; i < c.edgeCount * 2; ++i) {
      if (c.edgeProps.type[i] !== InteractionType.HydrogenBond) continue;
      ret.add(`${element(unit, c.a[i])}-${element(unit, c.b[i])}`);
    }
  }
  return ret;
}

describe('hydrogen bonds', async () => {
  // polar hydrogens (PDBQT)
  const parsed = await parsePDB(fs.readFileSync(path.join(examples, 'ace2.pdbqt'), 'utf8'), 'ace2', 'pdbqt').run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromPDB(parsed.result).run()).representative;
  const structure = Structure.ofModel(model);

  it('follow a donor hydrogen in another unit under another operator (per-atom symmetry)', async () => {
    const unit = structure.units[0] as Unit.Atomic;
    const { unitsContacts, unitsFeatures } = await computeInteractions(ctx, structure, {});
    const contacts = unitsContacts.get(unit.id);
    const { members, offsets, types } = unitsFeatures.get(unit.id);

    // an H-bond whose donor has exactly one hydrogen
    let donor = -1;
    let acceptor = -1;
    let hydrogen = -1;
    for (let i = 0; i < contacts.edgeCount * 2 && donor < 0; ++i) {
      if (contacts.edgeProps.type[i] !== InteractionType.HydrogenBond) continue;
      if (types[contacts.a[i]] !== FeatureType.HydrogenDonor) continue;
      const d = members[offsets[contacts.a[i]]] as StructureElement.UnitIndex;
      const hs: number[] = [];
      eachIntraBondedAtom(unit, d, (_, x) => {
        if (typeSymbol(unit, x) === 'H') hs.push(x);
      });
      if (hs.length !== 1) continue;
      donor = unit.elements[d];
      acceptor = unit.elements[members[offsets[contacts.b[i]]]];
      hydrogen = unit.elements[hs[0]];
    }
    expect(donor).toBeGreaterThanOrEqual(0);

    // the unit under a rotation R, whole and with that hydrogen in a unit of its own under R · (a rotation about the
    // hydrogen): the same position under another operator, bonded to the donor across units
    const R = Mat4.fromRotation(Mat4(), Math.PI / 5, Vec3.create(1, 2, 3));
    Mat4.setTranslation(R, Vec3.create(10, -20, 5));
    const xH = Vec3.create(
      model.atomicConformation.x[hydrogen],
      model.atomicConformation.y[hydrogen],
      model.atomicConformation.z[hydrogen],
    );
    const aboutH = Mat4.fromTranslation(Mat4(), xH);
    Mat4.mul(aboutH, aboutH, Mat4.fromRotation(Mat4(), Math.PI / 2, Vec3.unitZ));
    Mat4.mul(aboutH, aboutH, Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), xH)));
    const opR = SymmetryOperator.create('R', R);
    const opH = SymmetryOperator.create('R·H', Mat4.mul(Mat4(), R, aboutH));
    expect(Mat4.areEqual(opR.matrix, opH.matrix, 1e-3)).toBe(false);

    const create = (id: number, op: SymmetryOperator, elements: ArrayLike<number>) =>
      Unit.create(
        id,
        id,
        id,
        unit.traits,
        unit.kind,
        model,
        op,
        SortedArray.ofSortedArray<ElementIndex>(Int32Array.from(elements) as any),
      );
    const joined = Structure.create([create(0, opR, unit.elements)]);
    const split = Structure.create([
      create(
        0,
        opR,
        unit.elements.filter((e) => e !== hydrogen),
      ),
      create(1, opH, [hydrogen]),
    ]);
    expect(
      Vec3.distance(
        split.units[1].conformation.position(hydrogen as ElementIndex, Vec3()),
        joined.units[0].conformation.position(hydrogen as ElementIndex, Vec3()),
      ),
    ).toBeLessThan(1e-6);

    const expected = await hydrogenBonds(joined);
    expect(expected).toContain(`${donor}-${acceptor}`);
    expect(await hydrogenBonds(split)).toEqual(expected);
  });
});
