/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Fred Ludlow <Fred.Ludlow@astx.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { degToRad } from '../../../mol-math/misc';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { eachBondedAtom, typeSymbol } from './util';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';

/**
 * Numbering mostly inline with coordination number from VSEPR,
 * breaks with `SquarePlanar = 7`
 */
export const enum AtomGeometry {
    Spherical = 0,
    Terminal = 1,
    Linear = 2,
    Trigonal = 3,
    Tetrahedral = 4,
    TrigonalBiPyramidal = 5,
    Octahedral = 6,
    SquarePlanar = 7, // Okay, it breaks down somewhere!
    Unknown = 8
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

export function assignGeometry (totalCoordination: number): AtomGeometry {
    switch (totalCoordination) {
        case 0: return AtomGeometry.Spherical;
        case 1: return AtomGeometry.Terminal;
        case 2: return AtomGeometry.Linear;
        case 3: return AtomGeometry.Trigonal;
        case 4: return AtomGeometry.Tetrahedral;
        default: return AtomGeometry.Unknown;

    }
}

export const AtomGeometryAngles = new Map<AtomGeometry, number>([
    [ AtomGeometry.Linear, degToRad(180) ],
    [ AtomGeometry.Trigonal, degToRad(120) ],
    [ AtomGeometry.Tetrahedral, degToRad(109.4721) ],
    [ AtomGeometry.Octahedral, degToRad(90) ]
]);

// tmp objects for `calcAngles` and `calcPlaneAngle`
const tmpDir1 = Vec3();
const tmpDir2 = Vec3();
const tmpPosA = Vec3();
const tmpPosB = Vec3();
const tmpPosX = Vec3();

/**
 * Calculate the angles x-a1-a2 for all x where x is a heavy atom (not H) bonded to ap1.
 */
export function calcAngles (structure: Structure, unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex): number[] {
    const angles: number[] = [];
    unitA.conformation.position(unitA.elements[indexA], tmpPosA);
    unitB.conformation.position(unitB.elements[indexB], tmpPosB);
    Vec3.sub(tmpDir1, tmpPosB, tmpPosA);

    eachBondedAtom(structure, unitA, indexA, (unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) => {
        if (typeSymbol(unitX, indexX) !== Elements.H) {
            unitX.conformation.position(unitX.elements[indexX], tmpPosX);
            Vec3.sub(tmpDir2, tmpPosX, tmpPosA);
            angles.push(Vec3.angle(tmpDir1, tmpDir2));
        }
    });
    return angles;
}

/**
 * Find two neighbours of ap1 to define a plane (if possible) and
 * measure angle out of plane to ap2
 * @param  {AtomProxy} ap1 First atom (angle centre)
 * @param  {AtomProxy} ap2 Second atom (out-of-plane)
 * @return {number}        Angle from plane to second atom
 */
export function calcPlaneAngle (structure: Structure, unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex): number | undefined {
    unitA.conformation.position(unitA.elements[indexA], tmpPosA);
    unitB.conformation.position(unitB.elements[indexB], tmpPosB);
    Vec3.sub(tmpDir1, tmpPosB, tmpPosA);

    const neighbours = [Vec3(), Vec3()];
    let ni = 0;
    let unitX1: Unit.Atomic | undefined;
    let indexX1: StructureElement.UnitIndex | undefined;
    eachBondedAtom(structure, unitA, indexA, (unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) => {
        if (ni > 1) return;
        if (typeSymbol(unitX, indexX) !== Elements.H) {
            unitX1 = unitX;
            indexX1 = indexX;
            unitX.conformation.position(unitX.elements[indexX], tmpPosX);
            Vec3.sub(neighbours[ni++], tmpPosX, tmpPosA);
        }
    });
    if (ni === 1 && unitX1 && indexX1) {
        eachBondedAtom(structure, unitX1, indexX1, (unitX: Unit.Atomic, indexX: StructureElement.UnitIndex) => {
            if (ni > 1) return;
            if (unitX === unitA && indexX === indexA) return;
            if (typeSymbol(unitX, indexX) !== Elements.H) {
                unitX.conformation.position(unitX.elements[indexX], tmpPosX);
                Vec3.sub(neighbours[ni++], tmpPosX, tmpPosA);
            }
        });
    }

    if (ni !== 2) {
        return;
    }

    Vec3.cross(tmpDir2, neighbours[0], neighbours[1]);
    return Math.abs((Math.PI / 2) - Vec3.angle(tmpDir2, tmpDir1));
}