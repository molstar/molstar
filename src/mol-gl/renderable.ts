/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema, BaseValues } from './renderable/schema';
import { GraphicsRenderItem, ComputeRenderItem, GraphicsRenderVariant, MultiDrawInstancedData } from './webgl/render-item';
import { ValueCell } from '../mol-util';
import { idFactory } from '../mol-util/id-factory';
import { clamp } from '../mol-math/interpolate';
import { Frustum3D } from '../mol-math/geometry/primitives/frustum3d';
import { Plane3D } from '../mol-math/geometry/primitives/plane3d';
import { Sphere3D } from '../mol-math/geometry';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const p3distanceToPoint = Plane3D.distanceToPoint;
const f3intersectsSphere3D = Frustum3D.intersectsSphere3D;
const s3fromArray = Sphere3D.fromArray;

const getNextRenderableId = idFactory();

export type RenderableState = {
    disposed: boolean
    visible: boolean
    alphaFactor: number
    pickable: boolean
    colorOnly: boolean
    opaque: boolean
    writeDepth: boolean
}

export interface Renderable<T extends RenderableValues> {
    readonly id: number
    readonly materialId: number
    readonly values: T
    readonly state: RenderableState

    cull: (cameraPlane: Plane3D, frustum: Frustum3D) => void
    uncull: () => void
    render: (variant: GraphicsRenderVariant, sharedTexturesCount: number) => void
    getProgram: (variant: GraphicsRenderVariant) => Program
    update: () => void
    dispose: () => void
}

function getMdiData(cellCount: number, mdiData?: MultiDrawInstancedData): MultiDrawInstancedData {
    if (mdiData && mdiData.instanceCounts.length >= cellCount) {
        return mdiData;
    } else {
        return {
            firsts: new Int32Array(cellCount),
            counts: new Int32Array(cellCount),
            offsets: new Int32Array(cellCount),
            instanceCounts: new Int32Array(cellCount),
            baseVertices: new Int32Array(cellCount),
            baseInstances: new Uint32Array(cellCount),
            count: 0,
        };
    }
}

type GraphicsRenderableValues = RenderableValues & BaseValues

export function createRenderable<T extends GraphicsRenderableValues>(renderItem: GraphicsRenderItem, values: T, state: RenderableState): Renderable<T> {
    const id = getNextRenderableId();

    let mdiData = getMdiData(0);
    let cullEnabled = false;

    const s = Sphere3D();

    return {
        id,
        materialId: renderItem.materialId,
        values,
        state,

        cull: (cameraPlane: Plane3D, frustum: Frustum3D) => {
            cullEnabled = false;

            if (values.drawCount.ref.value === 0) return;
            if (values.instanceCount.ref.value === 0) return;
            if (!values.instanceGrid.ref.value) return;

            const { cellOffsets, cellSpheres, cellCount } = values.instanceGrid.ref.value;
            const [minDistance, maxDistance] = values.uLod.ref.value;
            const hasLod = minDistance !== 0 || maxDistance !== 0;

            mdiData = getMdiData(cellCount, mdiData);
            const { baseInstances, instanceCounts, counts } = mdiData;
            let o = 0;

            for (let i = 0; i < cellCount; ++i) {
                s3fromArray(s, cellSpheres, i * 4);
                if (hasLod) {
                    const d = p3distanceToPoint(cameraPlane, s.center);
                    if (d + s.radius < minDistance) continue;
                    if (d - s.radius > maxDistance) continue;
                }
                if (!f3intersectsSphere3D(frustum, s)) continue;

                const begin = cellOffsets[i];
                const end = cellOffsets[i + 1];
                const count = end - begin;
                if (count === 0) continue;

                if (o > 0 && baseInstances[o - 1] + instanceCounts[o - 1] === begin) {
                    instanceCounts[o - 1] += count;
                } else {
                    counts[o] = values.drawCount.ref.value;
                    instanceCounts[o] = count;
                    baseInstances[o] = begin;
                    o += 1;
                }
            }
            mdiData.count = o;

            // console.log({
            //     counts: counts.slice(),
            //     instanceCounts: instanceCounts.slice(),
            //     baseInstances: baseInstances.slice(),
            //     drawcount: mdiData.count,
            // });

            cullEnabled = true;
        },
        uncull: () => {
            cullEnabled = false;
        },
        render: (variant: GraphicsRenderVariant, sharedTexturesCount: number) => {
            if (values.uAlpha && values.alpha) {
                ValueCell.updateIfChanged(values.uAlpha, clamp(values.alpha.ref.value * state.alphaFactor, 0, 1));
            }
            renderItem.render(variant, sharedTexturesCount, cullEnabled ? mdiData : undefined);
        },
        getProgram: (variant: GraphicsRenderVariant) => renderItem.getProgram(variant),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    };
}

export type GraphicsRenderable = Renderable<GraphicsRenderableValues>

//

export interface ComputeRenderable<T extends RenderableValues> {
    readonly id: number
    readonly values: T

    render: () => void
    update: () => void
    dispose: () => void
}

export function createComputeRenderable<T extends Values<RenderableSchema>>(renderItem: ComputeRenderItem, values: T): ComputeRenderable<T> {
    return {
        id: getNextRenderableId(),
        values,

        render: () => renderItem.render('compute', 0),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    };
}