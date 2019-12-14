/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from '../../../mol-model/structure/structure';
import { ChunkedArray } from '../../../mol-data/util';
import { GridLookup3D } from '../../../mol-math/geometry';
import { OrderedSet } from '../../../mol-data/int';

export { Features }

interface Features {
    /** center x coordinate */
    readonly x: ArrayLike<number>
    /** center y coordinate */
    readonly y: ArrayLike<number>
    /** center z coordinate */
    readonly z: ArrayLike<number>
    /** number of features */
    readonly count: number
    readonly types: ArrayLike<FeatureType>
    readonly groups: ArrayLike<FeatureGroup>
    readonly offsets: ArrayLike<number>
    /** elements of this feature, range for feature i is offsets[i] to offsets[i + 1] */
    readonly members: ArrayLike<StructureElement.UnitIndex>
    /** lookup3d based on center coordinates */
    readonly lookup3d: GridLookup3D
}

namespace Features {
    /** maps elements to features, range for element i is offsets[i] to offsets[i + 1] */
    export type ElementsIndex = {
        readonly indices: ArrayLike<number>
        readonly offsets: ArrayLike<number>
    }

    export function createElementsIndex(features: Features, elementsCount: number): ElementsIndex {
        const offsets = new Int32Array(elementsCount + 1)
        const bucketFill = new Int32Array(elementsCount)
        const bucketSizes = new Int32Array(elementsCount)
        const { members, count, offsets: featureOffsets } = features
        for (let i = 0; i < count; ++i) ++bucketSizes[members[i]]

        let offset = 0
        for (let i = 0; i < elementsCount; i++) {
            offsets[i] = offset
            offset += bucketSizes[i]
        }
        offsets[elementsCount] = offset

        const indices = new Int32Array(offset)
        for (let i = 0; i < count; ++i) {
            for (let j = featureOffsets[i], jl = featureOffsets[i + 1]; j < jl; ++j) {
                const a = members[j]
                const oa = offsets[a] + bucketFill[a]
                indices[oa] = i
                ++bucketFill[a]
            }
        }

        return { indices, offsets }
    }
}

export { FeaturesBuilder }

interface FeaturesBuilder {
    clearState: () => void
    pushMember: (x: number, y: number, z: number, member: StructureElement.UnitIndex) => void
    addState: (type: FeatureType, group: FeatureGroup) => void
    addOne: (type: FeatureType, group: FeatureGroup, x: number, y: number, z: number, member: StructureElement.UnitIndex) => void
    getFeatures: () => Features
}

namespace FeaturesBuilder {
    interface State { x: number, y: number, z: number, offset: number, count: number }

    export function create(initialCount = 2048, chunkSize = 1024, features?: Features): FeaturesBuilder {
        const xCenters = ChunkedArray.create(Float32Array, 1, chunkSize, features ? features.x : initialCount)
        const yCenters = ChunkedArray.create(Float32Array, 1, chunkSize, features ? features.y : initialCount)
        const zCenters = ChunkedArray.create(Float32Array, 1, chunkSize, features ? features.z : initialCount)
        const types = ChunkedArray.create(Uint8Array, 1, chunkSize, features ? features.types : initialCount)
        const groups = ChunkedArray.create(Uint8Array, 1, chunkSize, features ? features.groups : initialCount)
        const offsets = ChunkedArray.create(Uint32Array, 1, chunkSize, features ? features.offsets : initialCount)
        const members = ChunkedArray.create(Uint32Array, 1, chunkSize, features ? features.members : initialCount)

        const state: State = { x: 0, y: 0, z: 0, offset: 0, count: 0 }

        return {
            clearState: () => {
                state.x = 0, state.y = 0, state.z = 0, state.offset = members.elementCount, state.count = 0
            },
            pushMember: (x: number, y: number, z: number, member: StructureElement.UnitIndex) => {
                ChunkedArray.add(members, member)
                state.x += x, state.y += y, state.z += z
            },
            addState: (type: FeatureType, group: FeatureGroup) => {
                const { count } = state
                if (count === 0) return
                ChunkedArray.add(types, type)
                ChunkedArray.add(groups, group)
                ChunkedArray.add(xCenters, state.x / count)
                ChunkedArray.add(yCenters, state.y / count)
                ChunkedArray.add(zCenters, state.z / count)
                ChunkedArray.add(offsets, state.offset)
            },
            addOne: (type: FeatureType, group: FeatureGroup, x: number, y: number, z: number, member: StructureElement.UnitIndex) => {
                ChunkedArray.add(types, type)
                ChunkedArray.add(groups, group)
                ChunkedArray.add(xCenters, x)
                ChunkedArray.add(yCenters, y)
                ChunkedArray.add(zCenters, z)
                ChunkedArray.add(offsets, members.elementCount)
                ChunkedArray.add(members, member)
            },
            getFeatures: () => {
                ChunkedArray.add(offsets, members.elementCount)
                const x = ChunkedArray.compact(xCenters, true) as ArrayLike<number>
                const y = ChunkedArray.compact(yCenters, true) as ArrayLike<number>
                const z = ChunkedArray.compact(zCenters, true) as ArrayLike<number>
                const count = xCenters.elementCount
                return {
                    x, y, z, count,
                    types: ChunkedArray.compact(types, true) as ArrayLike<FeatureType>,
                    groups: ChunkedArray.compact(groups, true) as ArrayLike<FeatureGroup>,
                    offsets: ChunkedArray.compact(offsets, true) as ArrayLike<number>,
                    members: ChunkedArray.compact(members, true) as ArrayLike<StructureElement.UnitIndex>,
                    lookup3d: GridLookup3D({ x, y, z, indices: OrderedSet.ofBounds(0, count) }),
                }
            }
        }
    }
}

export const enum FeatureType {
    None = 0,
    PositiveCharge = 1,
    NegativeCharge = 2,
    AromaticRing = 3,
    HydrogenDonor = 4,
    HydrogenAcceptor = 5,
    HalogenDonor = 6,
    HalogenAcceptor = 7,
    Hydrophobic = 8,
    WeakHydrogenDonor = 9,
    IonicTypePartner = 10,
    DativeBondPartner = 11,
    TransitionMetal = 12,
    IonicTypeMetal = 13
}

export const enum FeatureGroup {
    None = 0,
    QuaternaryAmine = 1,
    TertiaryAmine = 2,
    Sulfonium = 3,
    SulfonicAcid = 4,
    Sulfate = 5,
    Phosphate = 6,
    Halocarbon = 7,
    Guanidine = 8,
    Acetamidine = 9,
    Carboxylate = 10
}