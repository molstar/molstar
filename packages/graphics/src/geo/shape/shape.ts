/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Shape, ShapeGroup } from '@molstar/model/model/shape/shape';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry';
import { CentroidHelper } from '@molstar/core/math/geometry/centroid-helper';
import type { GroupMapping } from '../util.js';
import { ShapeGroupSizeTheme } from '@molstar/graphics/theme/size/shape-group';
import { ShapeGroupColorTheme } from '@molstar/graphics/theme/color/shape-group';
import type { Theme } from '@molstar/graphics/theme/theme';
import { type TransformData, createTransform as createGeometryTransform } from '../geometry/transform-data.js';
import { createRenderObject as createGlRenderObject, getNextMaterialId } from '@molstar/graphics/gl/render-object';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { LocationIterator } from '../util/location-iterator.js';
import { Geometry } from '../geometry/geometry.js';
import { OrderedSet } from '@molstar/core/data/int';

export function getTheme(shape: Shape) : Theme {
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
    return LocationIterator(shape.groupCount, instanceCount, 1, getLocation);
}

export function createTransform(transforms: Mat4[], invariantBoundingSphere: Sphere3D, cellSize: number, batchSize: number, transformData?: TransformData) {
    const transformArray = transformData && transformData.aTransform.ref.value.length >= transforms.length * 16 ? transformData.aTransform.ref.value : new Float32Array(transforms.length * 16);
    for (let i = 0, il = transforms.length; i < il; ++i) Mat4.toArray(transforms[i], transformArray, i * 16);
    return createGeometryTransform(transformArray, transforms.length, invariantBoundingSphere, cellSize, batchSize, transformData);
}

export function createRenderObject<G extends Geometry>(shape: Shape<G>, props: PD.Values<Geometry.Params<G>>) {
    const theme = getTheme(shape);
    const utils = Geometry.getUtils(shape.geometry);
    const materialId = getNextMaterialId();
    const locationIt = groupIterator(shape);
    const transform = createTransform(shape.transforms, shape.geometry.boundingSphere, props.cellSize, props.batchSize);
    const values = utils.createValues(shape.geometry, transform, locationIt, theme, props);
    const state = utils.createRenderableState(props);
    return createGlRenderObject(shape.geometry.kind, values, state, materialId);
}

const sphereHelper = new CentroidHelper(), tmpPos = Vec3.zero();

function sphereHelperInclude(groups: ShapeGroup.Loci['groups'], mapping: GroupMapping, positions: Float32Array, transforms: Mat4[]) {
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

function sphereHelperRadius(groups: ShapeGroup.Loci['groups'], mapping: GroupMapping, positions: Float32Array, transforms: Mat4[]) {
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

/** Graphics-side implementation for bounds of selected shape groups. */
export function getGroupBoundingSphere(loci: ShapeGroup.Loci, boundingSphere?: Sphere3D) {
    if (!boundingSphere) boundingSphere = Sphere3D();
    sphereHelper.reset();
    let padding = 0;
    const geometry = loci.shape.geometry as Geometry;
    const { transforms } = loci.shape;
    if (geometry.kind === 'mesh' || geometry.kind === 'points') {
        const positions = geometry.kind === 'mesh' ? geometry.vertexBuffer.ref.value : geometry.centerBuffer.ref.value;
        sphereHelperInclude(loci.groups, geometry.groupMapping, positions, transforms);
        sphereHelper.finishedIncludeStep();
        sphereHelperRadius(loci.groups, geometry.groupMapping, positions, transforms);
    } else if (geometry.kind === 'lines') {
        const start = geometry.startBuffer.ref.value, end = geometry.endBuffer.ref.value;
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
        for (const { ids, instance } of loci.groups) OrderedSet.forEach(ids, v => {
            const value = loci.shape.getSize(v, instance);
            if (padding < value) padding = value;
        });
    } else {
        return Sphere3D.copy(boundingSphere, geometry.boundingSphere);
    }
    Vec3.copy(boundingSphere.center, sphereHelper.center);
    boundingSphere.radius = Math.sqrt(sphereHelper.radiusSq);
    Sphere3D.expand(boundingSphere, boundingSphere, padding);
    return boundingSphere;
}

/** Create a shape with graphics-derived group counts and selected-group bounds. */
export function create<G extends Geometry>(name: string, sourceData: unknown, geometry: G, getColor: Shape['getColor'], getSize: Shape['getSize'], getLabel: Shape['getLabel'], transforms?: Mat4[], groupCount?: number): Shape<G> {
    const shape = Shape.create(name, sourceData, geometry, getColor, getSize, getLabel, groupCount ?? Geometry.getGroupCount(geometry), transforms);
    shape.getGroupBoundingSphere = (groups, out) => getGroupBoundingSphere(ShapeGroup.Loci(shape, groups), out);
    return shape;
}
