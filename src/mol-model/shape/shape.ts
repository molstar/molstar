/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { UUID } from '../../mol-util';
import { OrderedSet } from '../../mol-data/int';
import { Geometry } from '../../mol-geo/geometry/geometry';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Sphere3D } from '../../mol-math/geometry';
import { CentroidHelper } from '../../mol-math/geometry/centroid-helper';
import { GroupMapping } from '../../mol-geo/util';
import { ShapeGroupSizeTheme } from '../../mol-theme/size/shape-group';
import { ShapeGroupColorTheme } from '../../mol-theme/color/shape-group';
import { Theme } from '../../mol-theme/theme';
import { TransformData, createTransform as _createTransform } from '../../mol-geo/geometry/transform-data';
import { createRenderObject as _createRenderObject, getNextMaterialId } from '../../mol-gl/render-object';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { LocationIterator } from '../../mol-geo/util/location-iterator';

export interface Shape<G extends Geometry = Geometry> {
    /** A uuid to identify a shape object */
    readonly id: UUID
    /** A name to describe the shape */
    readonly name: string
    /** The data used to create the shape */
    readonly sourceData: unknown
    /** The geometry of the shape, e.g. `Mesh` or `Lines` */
    readonly geometry: G
    /** An array of transformation matrices to describe multiple instances of the geometry */
    readonly transforms: Mat4[]
    /** Number of groups in the geometry */
    readonly groupCount: number
    /** Get color for a given group */
    getColor(groupId: number, instanceId: number): Color
    /** Get size for a given group */
    getSize(groupId: number, instanceId: number): number
    /** Get label for a given group */
    getLabel(groupId: number, instanceId: number): string
}

export namespace Shape {
    export function create<G extends Geometry>(name: string, sourceData: unknown, geometry: G, getColor: Shape['getColor'], getSize: Shape['getSize'], getLabel: Shape['getLabel'], transforms?: Mat4[]): Shape<G> {
        return {
            id: UUID.create22(),
            name,
            sourceData,
            geometry,
            transforms: transforms || [Mat4.identity()],
            get groupCount() { return Geometry.getGroupCount(geometry); },
            getColor,
            getSize,
            getLabel
        };
    }

    export function getTheme(shape: Shape): Theme {
        return {
            color: ShapeGroupColorTheme({ shape }, {}),
            size: ShapeGroupSizeTheme({ shape }, {})
        };
    }

    export function groupIterator(shape: Shape): LocationIterator {
        const instanceCount = shape.transforms.length;
        const location = ShapeGroup.Location(shape);
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            location.group = groupIndex;
            location.instance = instanceIndex;
            return location;
        };
        return LocationIterator(shape.groupCount, instanceCount, getLocation);
    }

    export function createTransform(transforms: Mat4[], transformData?: TransformData) {
        const transformArray = transformData && transformData.aTransform.ref.value.length >= transforms.length * 16 ? transformData.aTransform.ref.value : new Float32Array(transforms.length * 16);
        for (let i = 0, il = transforms.length; i < il; ++i) {
            Mat4.toArray(transforms[i], transformArray, i * 16);
        }
        return _createTransform(transformArray, transforms.length, transformData);
    }

    export function createRenderObject<G extends Geometry>(shape: Shape<G>, props: PD.Values<Geometry.Params<G>>) {
        props;
        const theme = Shape.getTheme(shape);
        const utils = Geometry.getUtils(shape.geometry);

        const materialId = getNextMaterialId();
        const locationIt = groupIterator(shape);
        const transform = Shape.createTransform(shape.transforms);
        const values = utils.createValues(shape.geometry, transform, locationIt, theme, props);
        const state = utils.createRenderableState(props);

        return _createRenderObject(shape.geometry.kind, values, state, materialId);
    }

    export interface Loci { readonly kind: 'shape-loci', readonly shape: Shape }
    export function Loci(shape: Shape): Loci { return { kind: 'shape-loci', shape }; }
    export function isLoci(x: any): x is Loci { return !!x && x.kind === 'shape-loci'; }
    export function areLociEqual(a: Loci, b: Loci) { return a.shape === b.shape; }
    export function isLociEmpty(loci: Loci) { return loci.shape.groupCount === 0; }
}

export namespace ShapeGroup {
    export interface Location {
        readonly kind: 'group-location'
        shape: Shape
        group: number
        instance: number
    }

    export function Location(shape?: Shape, group = 0, instance = 0): Location {
        return { kind: 'group-location', shape: shape!, group, instance };
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'group-location';
    }

    export interface Loci {
        readonly kind: 'group-loci',
        readonly shape: Shape,
        readonly groups: ReadonlyArray<{
            readonly ids: OrderedSet<number>
            readonly instance: number
        }>
    }

    export function Loci(shape: Shape, groups: Loci['groups']): Loci {
        return { kind: 'group-loci', shape, groups: groups as Loci['groups'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'group-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.shape !== b.shape) return false;
        if (a.groups.length !== b.groups.length) return false;
        for (let i = 0, il = a.groups.length; i < il; ++i) {
            const { ids: idsA, instance: instanceA } = a.groups[i];
            const { ids: idsB, instance: instanceB } = b.groups[i];
            if (instanceA !== instanceB) return false;
            if (!OrderedSet.areEqual(idsA, idsB)) return false;
        }
        return true;
    }

    export function isLociEmpty(loci: Loci) {
        return size(loci) === 0 ? true : false;
    }

    export function size(loci: Loci) {
        let size = 0;
        for (const group of loci.groups) {
            size += OrderedSet.size(group.ids);
        }
        return size;
    }

    const sphereHelper = new CentroidHelper(), tmpPos = Vec3.zero();

    function sphereHelperInclude(groups: Loci['groups'], mapping: GroupMapping, positions: Float32Array, transforms: Mat4[]) {
        const { indices, offsets } = mapping;
        for (const { ids, instance } of groups) {
            OrderedSet.forEach(ids, v => {
                for (let i = offsets[v], il = offsets[v + 1]; i < il; ++i) {
                    Vec3.fromArray(tmpPos, positions, indices[i] * 3);
                    Vec3.transformMat4(tmpPos, tmpPos, transforms[instance]);
                    sphereHelper.includeStep(tmpPos);
                }
            });
        }
    }

    function sphereHelperRadius(groups: Loci['groups'], mapping: GroupMapping, positions: Float32Array, transforms: Mat4[]) {
        const { indices, offsets } = mapping;
        for (const { ids, instance } of groups) {
            OrderedSet.forEach(ids, v => {
                for (let i = offsets[v], il = offsets[v + 1]; i < il; ++i) {
                    Vec3.fromArray(tmpPos, positions, indices[i] * 3);
                    Vec3.transformMat4(tmpPos, tmpPos, transforms[instance]);
                    sphereHelper.radiusStep(tmpPos);
                }
            });
        }
    }

    export function getBoundingSphere(loci: Loci, boundingSphere?: Sphere3D) {
        if (!boundingSphere) boundingSphere = Sphere3D();

        sphereHelper.reset();
        let padding = 0;

        const { geometry, transforms } = loci.shape;

        if (geometry.kind === 'mesh' || geometry.kind === 'points') {
            const positions = geometry.kind === 'mesh'
                ? geometry.vertexBuffer.ref.value
                : geometry.centerBuffer.ref.value;
            sphereHelperInclude(loci.groups, geometry.groupMapping, positions, transforms);
            sphereHelper.finishedIncludeStep();
            sphereHelperRadius(loci.groups, geometry.groupMapping, positions, transforms);
        } else if (geometry.kind === 'lines') {
            const start = geometry.startBuffer.ref.value;
            const end = geometry.endBuffer.ref.value;
            sphereHelperInclude(loci.groups, geometry.groupMapping, start, transforms);
            sphereHelperInclude(loci.groups, geometry.groupMapping, end, transforms);
            sphereHelper.finishedIncludeStep();
            sphereHelperRadius(loci.groups, geometry.groupMapping, start, transforms);
            sphereHelperRadius(loci.groups, geometry.groupMapping, end, transforms);
        } else if (geometry.kind === 'spheres' || geometry.kind === 'text') {
            const positions = geometry.centerBuffer.ref.value;
            sphereHelperInclude(loci.groups, geometry.groupMapping, positions, transforms);
            sphereHelper.finishedIncludeStep();
            sphereHelperRadius(loci.groups, geometry.groupMapping, positions, transforms);
            for (const { ids, instance } of loci.groups) {
                OrderedSet.forEach(ids, v => {
                    const value = loci.shape.getSize(v, instance);
                    if (padding < value) padding = value;
                });
            }
        } else {
            // use whole shape bounding-sphere for other geometry kinds
            return Sphere3D.copy(boundingSphere, geometry.boundingSphere);
        }

        Vec3.copy(boundingSphere.center, sphereHelper.center);
        boundingSphere.radius = Math.sqrt(sphereHelper.radiusSq);
        Sphere3D.expand(boundingSphere, boundingSphere, padding);
        return boundingSphere;
    }
}