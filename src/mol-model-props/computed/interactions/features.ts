/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, Structure } from '../../../mol-model/structure/structure';
import { ChunkedArray } from '../../../mol-data/util';
import { GridLookup3D } from '../../../mol-math/geometry';
import { OrderedSet, SortedArray } from '../../../mol-data/int';
import { FeatureGroup, FeatureType } from './common';
import { ValenceModelProvider } from '../valence-model';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { getBoundary } from '../../../mol-math/geometry/boundary';

export { Features };

interface Features {
    /** number of features */
    readonly count: number
    /** center x coordinate, in invariant coordinate space */
    readonly x: ArrayLike<number>
    /** center y coordinate, in invariant coordinate space */
    readonly y: ArrayLike<number>
    /** center z coordinate, in invariant coordinate space */
    readonly z: ArrayLike<number>
    readonly types: ArrayLike<FeatureType>
    readonly groups: ArrayLike<FeatureGroup>
    readonly offsets: ArrayLike<number>
    /** elements of this feature, range for feature i is offsets[i] to offsets[i + 1] */
    readonly members: ArrayLike<StructureElement.UnitIndex>

    /** lookup3d based on center coordinates, in invariant coordinate space */
    readonly lookup3d: GridLookup3D<Features.FeatureIndex>
    /** maps unit elements to features, range for unit element i is offsets[i] to offsets[i + 1] */
    readonly elementsIndex: Features.ElementsIndex

    subset(types: ReadonlySet<FeatureType>): Features.Subset
}

namespace Features {
    /** Index into Features data arrays */
    export type FeatureIndex = { readonly '@type': 'feature-index' } & number

    export function setPosition(out: Vec3, unit: Unit, index: FeatureIndex, features: Features) {
        Vec3.set(out, features.x[index], features.y[index], features.z[index]);
        Vec3.transformMat4(out, out, unit.conformation.operator.matrix);
        return out;
    }

    /** maps unit elements to features, range for unit element i is offsets[i] to offsets[i + 1] */
    export type ElementsIndex = {
        /** feature indices */
        readonly indices: ArrayLike<FeatureIndex>
        /** range for unit element i is offsets[i] to offsets[i + 1] */
        readonly offsets: ArrayLike<number>
    }

    export type Data = {
        count: number
        x: ArrayLike<number>
        y: ArrayLike<number>
        z: ArrayLike<number>
        types: ArrayLike<FeatureType>
        groups: ArrayLike<FeatureGroup>
        offsets: ArrayLike<number>
        members: ArrayLike<StructureElement.UnitIndex>
    }

    export type Subset = {
        readonly indices: OrderedSet<FeatureIndex>
        readonly lookup3d: GridLookup3D
    }

    export function createElementsIndex(data: Data, elementsCount: number): ElementsIndex {
        const offsets = new Int32Array(elementsCount + 1);
        const bucketFill = new Int32Array(elementsCount);
        const bucketSizes = new Int32Array(elementsCount);
        const { members, count, offsets: featureOffsets } = data;
        for (let i = 0, il = featureOffsets[count]; i < il; ++i) ++bucketSizes[members[i]];

        let offset = 0;
        for (let i = 0; i < elementsCount; i++) {
            offsets[i] = offset;
            offset += bucketSizes[i];
        }
        offsets[elementsCount] = offset;

        const indices = new Int32Array(offset);
        for (let i = 0; i < count; ++i) {
            for (let j = featureOffsets[i], jl = featureOffsets[i + 1]; j < jl; ++j) {
                const a = members[j];
                const oa = offsets[a] + bucketFill[a];
                indices[oa] = i;
                ++bucketFill[a];
            }
        }

        return { indices: indices as unknown as ArrayLike<FeatureIndex>, offsets };
    }

    export function create(elementsCount: number, data: Data): Features {
        let lookup3d: GridLookup3D<FeatureIndex>;
        let elementsIndex: ElementsIndex;

        return {
            ...data,
            get lookup3d() {
                if (!lookup3d) {
                    const position = { x: data.x, y: data.y, z: data.z, indices: OrderedSet.ofBounds(0 as FeatureIndex, data.count as FeatureIndex) };
                    lookup3d = GridLookup3D(position, getBoundary(position));
                }
                return lookup3d;
            },
            get elementsIndex() {
                return elementsIndex || (elementsIndex = createElementsIndex(data, elementsCount));
            },

            subset: (types: Set<FeatureType>) => createSubset(data, types)
        };
    }

    export function createSubset(data: Data, types: ReadonlySet<FeatureType>): Subset {
        let lookup3d: GridLookup3D;

        const { count, types: _types } = data;
        const _indices = [];
        for (let i = 0; i < count; ++i) {
            if (types.has(_types[i])) _indices.push(i);
        }
        const indices = SortedArray.ofSortedArray<FeatureIndex>(_indices);

        return {
            indices,
            get lookup3d() {
                if (!lookup3d) {
                    const position = { x: data.x, y: data.y, z: data.z, indices };
                    lookup3d = GridLookup3D(position, getBoundary(position));
                }
                return lookup3d;
            }
        };
    }

    export interface Info {
        unit: Unit.Atomic,
        types: ArrayLike<FeatureType>,
        feature: FeatureIndex,
        x: ArrayLike<number>
        y: ArrayLike<number>
        z: ArrayLike<number>
        members: ArrayLike<StructureElement.UnitIndex>,
        offsets: ArrayLike<number>,
        idealGeometry: Int8Array
    }
    export function Info(structure: Structure, unit: Unit.Atomic, features: Features): Info {
        const valenceModel = ValenceModelProvider.get(structure).value;
        if (!valenceModel || !valenceModel.has(unit.id)) throw new Error('valence model required');

        return {
            unit,
            types: features.types,
            feature: -1 as any,
            x: features.x,
            y: features.y,
            z: features.z,
            members: features.members,
            offsets: features.offsets,
            idealGeometry: valenceModel.get(unit.id)!.idealGeometry
        };
    }

    export function position(out: Vec3, info: Info) {
        Vec3.set(out, info.x[info.feature], info.y[info.feature], info.z[info.feature]);
        Vec3.transformMat4(out, out, info.unit.conformation.operator.matrix);
        return out;
    }

    const tmpVecA = Vec3();
    const tmpVecB = Vec3();
    export function distance(infoA: Info, infoB: Info) {
        const elementA = infoA.members[infoA.offsets[infoA.feature]];
        const elementB = infoB.members[infoB.offsets[infoB.feature]];
        infoA.unit.conformation.position(infoA.unit.elements[elementA], tmpVecA);
        infoB.unit.conformation.position(infoB.unit.elements[elementB], tmpVecB);
        return Vec3.distance(tmpVecA, tmpVecB);
    }

    export interface Provider {
        types: Set<FeatureType>
        add: (structure: Structure, unit: Unit.Atomic, featuresBuilder: FeaturesBuilder) => void
    }
    export function Provider(types: FeatureType[], add: Provider['add']): Provider {
        return { types: new Set(types), add };
    }
}

export { FeaturesBuilder };

interface FeaturesBuilder {
    startState: () => void
    pushMember: (x: number, y: number, z: number, member: StructureElement.UnitIndex) => void
    finishState: (type: FeatureType, group: FeatureGroup) => void
    add: (type: FeatureType, group: FeatureGroup, x: number, y: number, z: number, member: StructureElement.UnitIndex) => void
    getFeatures: (elementsCount: number) => Features
}

namespace FeaturesBuilder {
    interface State { x: number, y: number, z: number, offset: number, count: number }

    export function create(initialCount = 2048, chunkSize = 1024, features?: Features): FeaturesBuilder {
        const xCenters = ChunkedArray.create(Float32Array, 1, chunkSize, features ? features.x : initialCount);
        const yCenters = ChunkedArray.create(Float32Array, 1, chunkSize, features ? features.y : initialCount);
        const zCenters = ChunkedArray.create(Float32Array, 1, chunkSize, features ? features.z : initialCount);
        const types = ChunkedArray.create(Uint8Array, 1, chunkSize, features ? features.types : initialCount);
        const groups = ChunkedArray.create(Uint8Array, 1, chunkSize, features ? features.groups : initialCount);
        const offsets = ChunkedArray.create(Uint32Array, 1, chunkSize, features ? features.offsets : initialCount);
        const members = ChunkedArray.create(Uint32Array, 1, chunkSize, features ? features.members : initialCount);

        const state: State = { x: 0, y: 0, z: 0, offset: 0, count: 0 };

        return {
            startState: () => {
                state.x = 0;
                state.y = 0;
                state.z = 0;
                state.offset = members.elementCount;
                state.count = 0;
            },
            pushMember: (x: number, y: number, z: number, member: StructureElement.UnitIndex) => {
                ChunkedArray.add(members, member);
                state.x += x;
                state.y += y;
                state.z += z;
                state.count += 1;
            },
            finishState: (type: FeatureType, group: FeatureGroup) => {
                const { count } = state;
                if (count === 0) return;
                ChunkedArray.add(types, type);
                ChunkedArray.add(groups, group);
                ChunkedArray.add(xCenters, state.x / count);
                ChunkedArray.add(yCenters, state.y / count);
                ChunkedArray.add(zCenters, state.z / count);
                ChunkedArray.add(offsets, state.offset);
            },
            add: (type: FeatureType, group: FeatureGroup, x: number, y: number, z: number, member: StructureElement.UnitIndex) => {
                ChunkedArray.add(types, type);
                ChunkedArray.add(groups, group);
                ChunkedArray.add(xCenters, x);
                ChunkedArray.add(yCenters, y);
                ChunkedArray.add(zCenters, z);
                ChunkedArray.add(offsets, members.elementCount);
                ChunkedArray.add(members, member);
            },
            getFeatures: (elementsCount: number) => {
                ChunkedArray.add(offsets, members.elementCount);
                const x = ChunkedArray.compact(xCenters, true) as ArrayLike<number>;
                const y = ChunkedArray.compact(yCenters, true) as ArrayLike<number>;
                const z = ChunkedArray.compact(zCenters, true) as ArrayLike<number>;
                const count = xCenters.elementCount;
                return Features.create(elementsCount, {
                    x, y, z, count,
                    types: ChunkedArray.compact(types, true) as ArrayLike<FeatureType>,
                    groups: ChunkedArray.compact(groups, true) as ArrayLike<FeatureGroup>,
                    offsets: ChunkedArray.compact(offsets, true) as ArrayLike<number>,
                    members: ChunkedArray.compact(members, true) as ArrayLike<StructureElement.UnitIndex>,
                });
            }
        };
    }
}