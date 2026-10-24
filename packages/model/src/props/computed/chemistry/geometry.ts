/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Fred Ludlow <Fred.Ludlow@astx.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { degToRad } from '@molstar/core/math/misc';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Structure, Unit, StructureElement } from '@molstar/model/model/structure';
import { eachBondedAtom, eachIntraBondedAtom, typeSymbol } from './util.js';
import { Elements } from '@molstar/model/model/structure/model/properties/atomic/types';

/**
 * Numbering mostly inline with coordination number from VSEPR,
 * breaks with `SquarePlanar = 7`
 */
export enum AtomGeometry {
  Spherical = 0,
  Terminal = 1,
  Linear = 2,
  Trigonal = 3,
  Tetrahedral = 4,
  TrigonalBiPyramidal = 5,
  Octahedral = 6,
  SquarePlanar = 7, // Okay, it breaks down somewhere!
  Unknown = 8,
}

export function geometryLabel(geometry: AtomGeometry): string {
  switch (geometry) {
    case AtomGeometry.Spherical:
      return 'Spherical';
    case AtomGeometry.Terminal:
      return 'Terminal';
    case AtomGeometry.Linear:
      return 'Linear';
    case AtomGeometry.Trigonal:
      return 'Trigonal';
    case AtomGeometry.Tetrahedral:
      return 'Tetrahedral';
    case AtomGeometry.TrigonalBiPyramidal:
      return 'Trigonal Bi-Pyramidal';
    case AtomGeometry.Octahedral:
      return 'Octahedral';
    case AtomGeometry.SquarePlanar:
      return 'Square Planar';
    case AtomGeometry.Unknown:
      return 'Unknown';
  }
}

export function assignGeometry(totalCoordination: number): AtomGeometry {
  switch (totalCoordination) {
    case 0:
      return AtomGeometry.Spherical;
    case 1:
      return AtomGeometry.Terminal;
    case 2:
      return AtomGeometry.Linear;
    case 3:
      return AtomGeometry.Trigonal;
    case 4:
      return AtomGeometry.Tetrahedral;
    default:
      return AtomGeometry.Unknown;
  }
}

export const AtomGeometryAngles = new Map<AtomGeometry, number>([
  [AtomGeometry.Linear, degToRad(180)],
  [AtomGeometry.Trigonal, degToRad(120)],
  [AtomGeometry.Tetrahedral, degToRad(109.4721)],
  [AtomGeometry.Octahedral, degToRad(90)],
]);

// tmp objects for `calcAngles` and `calcPlaneAngle`
const tmpDir1 = Vec3();
const tmpDir2 = Vec3();
const tmpPosA = Vec3();
const tmpPosB = Vec3();
const tmpPosX = Vec3();

/** Angles into the given arrays, through module state and callbacks: no closure per call (perf-audit-2 06/P12) */
const anglesState = {
  angles: [] as number[],
  anglesH: [] as number[],
  ignoreHydrogens: true,
};
function anglesVisitor(unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) {
  if (typeSymbol(unitX, indexX) !== Elements.H) {
    unitX.conformation.position(unitX.elements[indexX], tmpPosX);
    Vec3.sub(tmpDir2, tmpPosX, tmpPosA);
    anglesState.angles.push(Vec3.angle(tmpDir1, tmpDir2));
  } else if (!anglesState.ignoreHydrogens) {
    unitX.conformation.position(unitX.elements[indexX], tmpPosX);
    Vec3.sub(tmpDir2, tmpPosX, tmpPosA);
    anglesState.anglesH.push(Vec3.angle(tmpDir1, tmpDir2));
  }
}

/** `calcAngles` into `angles` and `anglesH` (emptied first) */
export function calcAnglesInto(
  structure: Structure,
  unitA: Unit.Atomic,
  indexA: StructureElement.UnitIndex,
  unitB: Unit.Atomic,
  indexB: StructureElement.UnitIndex,
  ignoreHydrogens: boolean,
  angles: number[],
  anglesH: number[],
) {
  angles.length = 0;
  anglesH.length = 0;
  unitA.conformation.position(unitA.elements[indexA], tmpPosA);
  unitB.conformation.position(unitB.elements[indexB], tmpPosB);
  Vec3.sub(tmpDir1, tmpPosB, tmpPosA);
  anglesState.angles = angles;
  anglesState.anglesH = anglesH;
  anglesState.ignoreHydrogens = ignoreHydrogens;
  eachBondedAtom(structure, unitA, indexA, anglesVisitor);
  anglesState.angles = anglesState.anglesH = NoAngles;
}
const NoAngles: number[] = [];

/**
 * Calculate the angles x-a1-a2 for all x where x is a heavy atom (not H) bonded to ap1.
 */
export function calcAngles(
  structure: Structure,
  unitA: Unit.Atomic,
  indexA: StructureElement.UnitIndex,
  unitB: Unit.Atomic,
  indexB: StructureElement.UnitIndex,
  ignoreHydrogens = true,
): [number[], number[]] {
  const angles: number[] = [];
  const anglesH: number[] = [];
  calcAnglesInto(structure, unitA, indexA, unitB, indexB, ignoreHydrogens, angles, anglesH);
  return [angles, anglesH];
}

const planeState = {
  neighbours: [Vec3(), Vec3()],
  ni: 0,
  unitA: undefined as Unit.Atomic | undefined,
  indexA: -1 as StructureElement.UnitIndex,
  unitX1: undefined as Unit.Atomic | undefined,
  indexX1: undefined as StructureElement.UnitIndex | undefined,
};
function planeFirstVisitor(unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) {
  const p = planeState;
  if (p.ni > 1) return;
  if (typeSymbol(unitX, indexX) !== Elements.H) {
    p.unitX1 = unitX;
    p.indexX1 = indexX;
    unitX.conformation.position(unitX.elements[indexX], tmpPosX);
    Vec3.sub(p.neighbours[p.ni++], tmpPosX, tmpPosA);
  }
}
function planeSecondVisitor(unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) {
  const p = planeState;
  if (p.ni > 1) return;
  if (unitX === p.unitA && indexX === p.indexA) return;
  if (typeSymbol(unitX, indexX) !== Elements.H) {
    unitX.conformation.position(unitX.elements[indexX], tmpPosX);
    Vec3.sub(p.neighbours[p.ni++], tmpPosX, tmpPosA);
  }
}

/**
 * Find two neighbours of ap1 to define a plane (if possible) and
 * measure angle out of plane to ap2
 * @param  {AtomProxy} ap1 First atom (angle centre)
 * @param  {AtomProxy} ap2 Second atom (out-of-plane)
 * @return {number}        Angle from plane to second atom
 */
export function calcPlaneAngle(
  structure: Structure,
  unitA: Unit.Atomic,
  indexA: StructureElement.UnitIndex,
  unitB: Unit.Atomic,
  indexB: StructureElement.UnitIndex,
): number | undefined {
  unitA.conformation.position(unitA.elements[indexA], tmpPosA);
  unitB.conformation.position(unitB.elements[indexB], tmpPosB);
  Vec3.sub(tmpDir1, tmpPosB, tmpPosA);

  const p = planeState;
  p.ni = 0;
  p.unitA = unitA;
  p.indexA = indexA;
  p.unitX1 = undefined;
  p.indexX1 = undefined;
  eachBondedAtom(structure, unitA, indexA, planeFirstVisitor);
  // index 0 is a neighbour too (perf-audit-2 06/B4)
  if (p.ni === 1 && p.unitX1 && p.indexX1 !== undefined) {
    eachBondedAtom(structure, p.unitX1, p.indexX1, planeSecondVisitor);
  }
  p.unitA = p.unitX1 = undefined;

  if (p.ni !== 2) {
    return;
  }

  Vec3.cross(tmpDir2, p.neighbours[0], p.neighbours[1]);
  return Math.abs(Math.PI / 2 - Vec3.angle(tmpDir2, tmpDir1));
}

const hydrogenState = {
  out: undefined as unknown as { unit: Unit.Atomic; index: StructureElement.UnitIndex },
  minDistSq: 0,
};
function closestHydrogenVisitor(unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) {
  if (typeSymbol(unitX, indexX) === Elements.H) {
    unitX.conformation.position(unitX.elements[indexX], tmpPosX);
    const dist = Vec3.squaredDistance(tmpPosX, tmpPosB);
    if (dist < hydrogenState.minDistSq) {
      hydrogenState.minDistSq = dist;
      hydrogenState.out.unit = unitX;
      hydrogenState.out.index = indexX;
    }
  }
}

/** `closestHydrogen` into `out` */
export function closestHydrogenInto(
  structure: Structure,
  unitA: Unit.Atomic,
  indexA: StructureElement.UnitIndex,
  unitB: Unit.Atomic,
  indexB: StructureElement.UnitIndex,
  out: { unit: Unit.Atomic; index: StructureElement.UnitIndex },
) {
  out.unit = unitA;
  out.index = indexA;
  unitA.conformation.position(unitA.elements[indexA], tmpPosA);
  unitB.conformation.position(unitB.elements[indexB], tmpPosB);
  hydrogenState.out = out;
  hydrogenState.minDistSq = Vec3.squaredDistance(tmpPosA, tmpPosB);
  eachBondedAtom(structure, unitA, indexA, closestHydrogenVisitor);
  // the module state keeps no unit (and its model) after the call
  hydrogenState.out = undefined as unknown as typeof out;
  return out;
}

/**
 * The hydrogen covalently bonded to (unitA, indexA) nearest to (unitB, indexB), or (unitA, indexA) itself if none is
 * nearer. Follows bonds across units: with per-atom symmetry (small-molecule crystals) a hydrogen can be in another
 * unit, under another operator, than the atom it is bonded to.
 */
export function closestHydrogen(
  structure: Structure,
  unitA: Unit.Atomic,
  indexA: StructureElement.UnitIndex,
  unitB: Unit.Atomic,
  indexB: StructureElement.UnitIndex,
): { unit: Unit.Atomic; index: StructureElement.UnitIndex } {
  return closestHydrogenInto(structure, unitA, indexA, unitB, indexB, { unit: unitA, index: indexA });
}

/**
 * The index in `unitA` of the hydrogen bonded to (unitA, indexA) nearest to (unitB, indexB), or `indexA` if none is
 * nearer. Only hydrogens of `unitA`, so the index is always one of `unitA`'s; `closestHydrogen` also follows bonds
 * across units.
 */
export function closestHydrogenIndex(
  _structure: Structure,
  unitA: Unit.Atomic,
  indexA: StructureElement.UnitIndex,
  unitB: Unit.Atomic,
  indexB: StructureElement.UnitIndex,
) {
  let hIndex = indexA;
  unitA.conformation.position(unitA.elements[indexA], tmpPosA);
  unitB.conformation.position(unitB.elements[indexB], tmpPosB);
  Vec3.sub(tmpDir1, tmpPosB, tmpPosA);
  let minDistSq = Vec3.squaredDistance(tmpPosA, tmpPosB);
  eachIntraBondedAtom(unitA, indexA, (unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) => {
    if (typeSymbol(unitX, indexX) === Elements.H) {
      unitX.conformation.position(unitX.elements[indexX], tmpPosX);
      const dist = Vec3.squaredDistance(tmpPosX, tmpPosB);
      if (dist < minDistSq) {
        minDistSq = dist;
        hIndex = indexX;
      }
    }
  });
  return hIndex;
}
