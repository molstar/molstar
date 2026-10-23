/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SortedArray } from '@molstar/core/data/int';
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { Structure, type StructureElement, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { calcPlaneAngle, closestHydrogen, closestHydrogenIndex } from '../geometry.js';
import { eachBondedAtom, typeSymbol } from '../util.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

describe('closestHydrogen, closestHydrogenIndex', async () => {
  // polar hydrogens (PDBQT)
  const parsed = await parsePDB(fs.readFileSync(path.join(examples, 'ace2.pdbqt'), 'utf8'), 'ace2', 'pdbqt').run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromPDB(parsed.result).run()).representative;
  const structure = Structure.ofModel(model);
  const unit = structure.units[0] as Unit.Atomic;

  // the first atom with a bonded hydrogen, and that hydrogen
  let donor = -1 as StructureElement.UnitIndex;
  let hydrogen = -1 as StructureElement.UnitIndex;
  for (let i = 0; i < unit.elements.length && donor < 0; ++i) {
    const index = i as StructureElement.UnitIndex;
    if (typeSymbol(unit, index) === 'H') continue;
    eachBondedAtom(structure, unit, index, (_, x) => {
      if (donor < 0 && typeSymbol(unit, x) === 'H') {
        donor = index;
        hydrogen = x;
      }
    });
  }

  it("returns the donor unit's hydrogen nearest the other atom", () => {
    expect(donor).toBeGreaterThanOrEqual(0);
    // the other atom is the hydrogen itself: nearest at distance 0
    expect(closestHydrogenIndex(structure, unit, donor, unit, hydrogen)).toBe(hydrogen);
  });

  it("closestHydrogenIndex ignores a hydrogen bonded across units (the index is always the donor unit's), closestHydrogen follows it", () => {
    // the same hydrogen in a unit of its own, covalently bonded to the donor across units
    const e = unit.elements[hydrogen];
    const rest = SortedArray.ofSortedArray<ElementIndex>(Int32Array.from(unit.elements.filter((x) => x !== e)) as any);
    const heavy = Unit.create(0, 0, 0, unit.traits, unit.kind, unit.model, unit.conformation.operator, rest);
    const single = SortedArray.ofSingleton(e) as SortedArray<ElementIndex>;
    const h = Unit.create(
      1,
      1,
      1,
      unit.traits,
      unit.kind,
      unit.model,
      unit.conformation.operator,
      single,
    ) as Unit.Atomic;
    const split = Structure.create([heavy, h]);
    const donorInHeavy = SortedArray.indexOf(heavy.elements, unit.elements[donor]) as StructureElement.UnitIndex;

    let bondedAcross = false;
    eachBondedAtom(split, heavy as Unit.Atomic, donorInHeavy, (u) => {
      if (u === h) bondedAcross = true;
    });
    expect(bondedAcross).toBe(true);

    const index = closestHydrogenIndex(split, heavy as Unit.Atomic, donorInHeavy, h, 0 as StructureElement.UnitIndex);
    expect(index).toBe(donorInHeavy);
    // closestHydrogen follows the bond across units and says which unit the hydrogen is in
    const closest = closestHydrogen(split, heavy as Unit.Atomic, donorInHeavy, h, 0 as StructureElement.UnitIndex);
    expect(closest.unit).toBe(h);
    expect(closest.index).toBe(0);
  });
});

describe('calcPlaneAngle', () => {
  it('finds the plane through a neighbour at unit index 0 (perf-audit-2 06/B4)', async () => {
    // acetate with its carbon first, and a water above the carboxyl plane
    const pdb = [
      'HETATM    1  C   ACT A   1       0.000   0.000   0.000  1.00  0.00           C',
      'HETATM    2  O   ACT A   1       1.250   0.000   0.000  1.00  0.00           O',
      'HETATM    3  OXT ACT A   1      -0.625   1.083   0.000  1.00  0.00           O',
      'HETATM    4  CH3 ACT A   1      -0.760  -1.316   0.000  1.00  0.00           C',
      'HETATM    5  O   HOH B   2       3.000   0.000   1.000  1.00  0.00           O',
      'END',
    ].join('\n');
    const parsed = await parsePDB(pdb, 'act').run();
    if (parsed.isError) throw new Error(parsed.message);
    const structure = Structure.ofModel((await trajectoryFromPDB(parsed.result).run()).representative);
    const act = structure.units.find((u) => u.elements.length === 4) as Unit.Atomic;
    const water = structure.units.find((u) => u.elements.length === 1) as Unit.Atomic;
    expect(typeSymbol(act, 0 as StructureElement.UnitIndex)).toBe('C');

    // O's only heavy neighbour is the carbon at index 0; the plane takes the carbon's other neighbour
    const angle = calcPlaneAngle(
      structure,
      act,
      1 as StructureElement.UnitIndex,
      water,
      0 as StructureElement.UnitIndex,
    );
    expect(angle).toBeDefined();
    expect(angle!).toBeGreaterThan(0.1);
  });
});
