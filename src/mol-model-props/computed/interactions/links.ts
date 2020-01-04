/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit } from '../../../mol-model/structure';
import { Features } from './features';
import { InteractionType, FeatureType } from './common';
import { IntraLinksBuilder, InterLinksBuilder } from './links-builder';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { altLoc, connectedTo } from '../chemistry/util';
import { OrderedSet } from '../../../mol-data/int';

const MAX_DISTANCE = 5

export interface LinkProvider<P extends PD.Params> {
    readonly name: string
    readonly params: P
    createTester(props: PD.Values<P>): LinkTester
}

export interface LinkTester {
    readonly maxDistance: number
    readonly requiredFeatures: ReadonlySet<FeatureType>
    getType: (structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number) => InteractionType | undefined
}

function validPair(structure: Structure, infoA: Features.Info, infoB: Features.Info): boolean {
    const indexA = infoA.members[infoA.offsets[infoA.feature]]
    const indexB = infoB.members[infoB.offsets[infoB.feature]]
    if (indexA === indexB) return false // no self interaction
    const altA = altLoc(infoA.unit, indexA)
    const altB = altLoc(infoB.unit, indexB)
    if (altA && altB && altA !== altB) return false // incompatible alternate location id
    if (infoA.unit.residueIndex[infoA.unit.elements[indexA]] === infoB.unit.residueIndex[infoB.unit.elements[indexB]]) return false // same residue
    // no hbond if donor and acceptor are bonded
    if (connectedTo(structure, infoA.unit, indexA, infoB.unit, indexB)) return false

    return true
}

/**
 * Add all intra-unit links, i.e. pairs of features
 */
export function addUnitLinks(structure: Structure, unit: Unit.Atomic, features: Features, builder: IntraLinksBuilder, testers: ReadonlyArray<LinkTester>) {

    for (const tester of testers) {
        _addUnitLinks(structure, unit, features, builder, tester)
    }
}

function _addUnitLinks(structure: Structure, unit: Unit.Atomic, features: Features, builder: IntraLinksBuilder, tester: LinkTester) {
    const { x, y, z } = features
    const { lookup3d, indices: subsetIndices } = features.subset(tester.requiredFeatures)

    const infoA = Features.Info(structure, unit, features)
    const infoB = { ...infoA }

    for (let t = 0, tl = OrderedSet.size(subsetIndices); t < tl; ++t) {
        const i = OrderedSet.getAt(subsetIndices, t)
        const { count, indices, squaredDistances } = lookup3d.find(x[i], y[i], z[i], tester.maxDistance)
        if (count === 0) continue

        infoA.feature = i

        for (let r = 0; r < count; ++r) {
            const j = OrderedSet.getAt(subsetIndices, indices[r])
            if (j <= i) continue

            infoB.feature = j
            if (!validPair(structure, infoA, infoB)) continue

            const type = tester.getType(structure, infoA, infoB, squaredDistances[r])
            if (type) builder.add(i, j, type)
        }
    }
}

const _imageTransform = Mat4()

/**
 * Add all inter-unit links, i.e. pairs of features
 */
export function addStructureLinks(structure: Structure, unitA: Unit.Atomic, featuresA: Features, unitB: Unit.Atomic, featuresB: Features, builder: InterLinksBuilder, testers: ReadonlyArray<LinkTester>) {
    const { count, x: xA, y: yA, z: zA } = featuresA;
    const { lookup3d } = featuresB;

    // the lookup queries need to happen in the "unitB space".
    // that means imageA = inverseOperB(operA(i))
    const imageTransform = Mat4.mul(_imageTransform, unitB.conformation.operator.inverse, unitA.conformation.operator.matrix)
    const isNotIdentity = !Mat4.isIdentity(imageTransform)
    const imageA = Vec3()

    const maxDistance = Math.max(...testers.map(t => t.maxDistance))
    const { center, radius } = lookup3d.boundary.sphere;
    const testDistanceSq = (radius + maxDistance) * (radius + maxDistance);

    const infoA = Features.Info(structure, unitA, featuresA)
    const infoB = Features.Info(structure, unitB, featuresB)

    builder.startUnitPair(unitA, unitB)

    for (let i = 0; i < count; ++i) {
        Vec3.set(imageA, xA[i], yA[i], zA[i])
        if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform)
        if (Vec3.squaredDistance(imageA, center) > testDistanceSq) continue

        const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], MAX_DISTANCE)
        if (count === 0) continue

        infoA.feature = i

        for (let r = 0; r < count; ++r) {
            const j = indices[r]
            infoB.feature = j
            if (!validPair(structure, infoA, infoB)) continue

            const distanceSq = squaredDistances[r]
            for (const tester of testers) {
                if (distanceSq < tester.maxDistance * tester.maxDistance) {
                    const type = tester.getType(structure, infoA, infoB, distanceSq)
                    if (type) {
                        builder.add(i, j, type)
                        break
                    }
                }
            }
        }
    }

    builder.finishUnitPair()
}