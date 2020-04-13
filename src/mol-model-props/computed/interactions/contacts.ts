/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { Features } from './features';
import { InteractionType, FeatureType } from './common';
import { IntraContactsBuilder, InterContactsBuilder } from './contacts-builder';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { altLoc, connectedTo, typeSymbol } from '../chemistry/util';
import { OrderedSet } from '../../../mol-data/int';
import { VdwRadius } from '../../../mol-model/structure/model/properties/atomic';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';

export const ContactsParams = {
    lineOfSightDistFactor: PD.Numeric(1.0, { min: 0, max: 3, step: 0.1 }),
};
export type ContactsParams = typeof ContactsParams
export type ContactsProps = PD.Values<ContactsParams>

const MAX_LINE_OF_SIGHT_DISTANCE = 3;

export interface ContactProvider<P extends PD.Params> {
    readonly name: string
    readonly params: P
    createTester(props: PD.Values<P>): ContactTester
}

export interface ContactTester {
    readonly maxDistance: number
    readonly requiredFeatures: ReadonlySet<FeatureType>
    getType: (structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number) => InteractionType | undefined
}

function validPair(structure: Structure, infoA: Features.Info, infoB: Features.Info): boolean {
    const indexA = infoA.members[infoA.offsets[infoA.feature]];
    const indexB = infoB.members[infoB.offsets[infoB.feature]];
    if (indexA === indexB) return false; // no self interaction
    const altA = altLoc(infoA.unit, indexA);
    const altB = altLoc(infoB.unit, indexB);
    if (altA && altB && altA !== altB) return false; // incompatible alternate location id
    if (infoA.unit.residueIndex[infoA.unit.elements[indexA]] === infoB.unit.residueIndex[infoB.unit.elements[indexB]]) return false; // same residue
    // e.g. no hbond if donor and acceptor are bonded
    if (connectedTo(structure, infoA.unit, indexA, infoB.unit, indexB)) return false;

    return true;
}

//

function invalidAltLoc (unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    const altA = altLoc(unitA, indexA);
    const altB = altLoc(unitB, indexB);
    return altA && altB && altA !== altB;
}

function isMember(element: StructureElement.UnitIndex, info: Features.Info) {
    const { feature, offsets, members } = info;
    for (let i = offsets[feature], il = offsets[feature + 1]; i < il; ++i) {
        if (members[i] === element) return true;
    }
    return false;
}

const tmpVec = Vec3();
const tmpVecA = Vec3();
const tmpVecB = Vec3();

function checkLineOfSight(structure: Structure, infoA: Features.Info, infoB: Features.Info, distFactor: number) {
    const featureA = infoA.feature;
    const featureB = infoB.feature;
    const indexA = infoA.members[infoA.offsets[featureA]];
    const indexB = infoB.members[infoB.offsets[featureB]];

    Features.position(tmpVecA, infoA);
    Features.position(tmpVecB, infoB);
    Vec3.scale(tmpVec, Vec3.add(tmpVec, tmpVecA, tmpVecB), 0.5);

    const distMax = distFactor * MAX_LINE_OF_SIGHT_DISTANCE;

    const { count, indices, units, squaredDistances } = structure.lookup3d.find(tmpVec[0], tmpVec[1], tmpVec[2], distMax);
    if (count === 0) return true;

    for (let r = 0; r < count; ++r) {
        const i = indices[r];
        const unit = units[r];
        if (!Unit.isAtomic(unit)) continue;

        const element = typeSymbol(unit, i);
        // allow hydrogens
        if (element === Elements.H) continue;

        const vdw = VdwRadius(element);
        // check distance
        if (vdw * vdw * distFactor * distFactor <= squaredDistances[r]) continue;

        // allow different altlocs
        if (invalidAltLoc(unit, i, infoA.unit, indexA) || invalidAltLoc(unit, i, infoB.unit, indexB)) continue;

        // allow member atoms
        if ((infoA.unit === unit && isMember(i, infoA)) || (infoB.unit === unit && isMember(i, infoB))) continue;

        unit.conformation.position(unit.elements[i], tmpVec);
        // allow atoms at the center of functional groups
        if (Vec3.squaredDistance(tmpVec, tmpVecA) < 1 || Vec3.squaredDistance(tmpVec, tmpVecB) < 1) continue;

        return false;
    }

    return true;
}

/**
 * Add all intra-unit contacts, i.e. pairs of features
 */
export function addUnitContacts(structure: Structure, unit: Unit.Atomic, features: Features, builder: IntraContactsBuilder, testers: ReadonlyArray<ContactTester>, props: ContactsProps) {
    for (const tester of testers) {
        _addUnitContacts(structure, unit, features, builder, tester, props);
    }
}

function _addUnitContacts(structure: Structure, unit: Unit.Atomic, features: Features, builder: IntraContactsBuilder, tester: ContactTester, props: ContactsProps) {
    const { x, y, z } = features;
    const { lookup3d, indices: subsetIndices } = features.subset(tester.requiredFeatures);

    const infoA = Features.Info(structure, unit, features);
    const infoB = { ...infoA };

    const distFactor = props.lineOfSightDistFactor;

    for (let t = 0, tl = OrderedSet.size(subsetIndices); t < tl; ++t) {
        const i = OrderedSet.getAt(subsetIndices, t);
        const { count, indices, squaredDistances } = lookup3d.find(x[i], y[i], z[i], tester.maxDistance);
        if (count === 0) continue;

        infoA.feature = i;

        for (let r = 0; r < count; ++r) {
            const j = OrderedSet.getAt(subsetIndices, indices[r]);
            if (j <= i) continue;

            infoB.feature = j;
            if (!validPair(structure, infoA, infoB)) continue;

            const type = tester.getType(structure, infoA, infoB, squaredDistances[r]);
            if (type && checkLineOfSight(structure, infoA, infoB, distFactor)) {
                builder.add(i, j, type);
            }
        }
    }
}

const _imageTransform = Mat4();

/**
 * Add all inter-unit contacts, i.e. pairs of features
 */
export function addStructureContacts(structure: Structure, unitA: Unit.Atomic, featuresA: Features, unitB: Unit.Atomic, featuresB: Features, builder: InterContactsBuilder, testers: ReadonlyArray<ContactTester>, props: ContactsProps) {
    const { count: countA, x: xA, y: yA, z: zA } = featuresA;
    const { lookup3d } = featuresB;

    // the lookup queries need to happen in the "unitB space".
    // that means imageA = inverseOperB(operA(i))
    const imageTransform = Mat4.mul(_imageTransform, unitB.conformation.operator.inverse, unitA.conformation.operator.matrix);
    const isNotIdentity = !Mat4.isIdentity(imageTransform);
    const imageA = Vec3();

    const maxDistance = Math.max(...testers.map(t => t.maxDistance));
    const { center, radius } = lookup3d.boundary.sphere;
    const testDistanceSq = (radius + maxDistance) * (radius + maxDistance);

    const distFactor = props.lineOfSightDistFactor;

    const infoA = Features.Info(structure, unitA, featuresA);
    const infoB = Features.Info(structure, unitB, featuresB);

    builder.startUnitPair(unitA, unitB);

    for (let i = 0 as Features.FeatureIndex; i < countA; ++i) {
        Vec3.set(imageA, xA[i], yA[i], zA[i]);
        if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform);
        if (Vec3.squaredDistance(imageA, center) > testDistanceSq) continue;

        const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], maxDistance);
        if (count === 0) continue;

        infoA.feature = i;

        for (let r = 0; r < count; ++r) {
            const j = indices[r];
            infoB.feature = j;
            if (!validPair(structure, infoA, infoB)) continue;

            const distanceSq = squaredDistances[r];
            for (const tester of testers) {
                if (distanceSq < tester.maxDistance * tester.maxDistance) {
                    const type = tester.getType(structure, infoA, infoB, distanceSq);
                    if (type && checkLineOfSight(structure, infoA, infoB, distFactor)) {
                        builder.add(i, j, type);
                        break;
                    }
                }
            }
        }
    }

    builder.finishUnitPair();
}