/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema, BaseValues } from './renderable/schema';
import { GraphicsRenderItem, ComputeRenderItem, GraphicsRenderVariant, MultiDrawBaseData, Transparency } from './webgl/render-item';
import { ValueCell } from '../mol-util/value-cell';
import { idFactory } from '../mol-util/id-factory';
import { clamp } from '../mol-math/interpolate';
import { Frustum3D } from '../mol-math/geometry/primitives/frustum3d';
import { Plane3D } from '../mol-math/geometry/primitives/plane3d';
import { Sphere3D } from '../mol-math/geometry/primitives/sphere3d';
import { Vec4 } from '../mol-math/linear-algebra/3d/vec4';
import { WebGLStats } from './webgl/context';
import { isTimingMode } from '../mol-util/debug';
import { InstanceGrid } from '../mol-math/geometry/instance-grid';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const p3distanceToPoint = Plane3D.distanceToPoint;
const f3intersectsSphere3D = Frustum3D.intersectsSphere3D;
const f3containsSphere3D = Frustum3D.containsSphere3D;
const s3fromArray = Sphere3D.fromArray;

const getNextRenderableId = idFactory();

// Opaque, unforgeable token identifying "a frame" for cull-cache reuse. Deliberately an
// object identity rather than a plain number: several independent owners (Canvas3D's main
// render loop, ImagePass, PickHelper, RayHelper, ...) each mint their own frame tokens while
// culling the SAME underlying renderables - if `frame` were a plain counter, two unrelated
// owners could coincidentally reach the same integer at the same time and wrongly reuse each
// other's cull result. Object identity (`===`) makes that structurally impossible: a token
// from one owner can never equal a token from another, no coordination between owners needed.
declare const frameBrand: unique symbol;
export type Frame = { readonly [frameBrand]: never };
/** Mint a new frame token - call this whenever an owner starts a new logical "frame". */
export function createFrame(): Frame {
    return {} as Frame;
}

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

    cull: (cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null, stats: WebGLStats, frame: Frame) => void
    uncull: () => void
    cullSimple: (d: number, radius: number, scale: number) => void
    render: (variant: GraphicsRenderVariant, sharedTexturesCount: number) => void
    getByteCount: () => number
    getProgram: (variant: GraphicsRenderVariant) => Program
    setTransparency: (transparency: Transparency) => void
    update: () => void
    dispose: () => void
}

function getMdbData(cellCount: number, mdbData?: MultiDrawBaseData): MultiDrawBaseData {
    if (mdbData && mdbData.instanceCounts.length >= cellCount) {
        return mdbData;
    } else {
        return {
            firsts: new Int32Array(cellCount),
            counts: new Int32Array(cellCount),
            offsets: new Int32Array(cellCount),
            instanceCounts: new Int32Array(cellCount),
            baseVertices: new Int32Array(cellCount),
            baseInstances: new Uint32Array(cellCount),
            count: 0,
            uniforms: [],
        };
    }
}

//

type LodLevelsValue = [minDistance: number, maxDistance: number, overlap: number, count: number, scale: number][]

/** The minimal set of values needed to cull a render item. */
type CullValues = {
    readonly drawCount: ValueCell<number>
    readonly instanceCount: ValueCell<number>
    readonly uLod: ValueCell<Vec4>
    readonly instanceGrid: ValueCell<InstanceGrid>
    readonly boundingSphere: ValueCell<Sphere3D>
    readonly lodLevels?: ValueCell<unknown>
}

// `frame` is bumped once per Canvas3D render() call (see canvas3d.ts), spanning both
// stereo eyes and all multi-sample jitter sub-renders of that call - a renderable's
// `update()` (the only thing that can change `values` mid-session) is only ever
// called before a frame's first `cull()`, never in between, so reusing the last cull
// result for the rest of a frame is exact, not an approximation; reusing ACROSS eyes/
// jitter samples within that same frame is a deliberate, accepted approximation.

/** A multi-draw-base-data list, one entry per lod level. */
interface MdbList {
    readonly list: MultiDrawBaseData[]
    /** resize and reset the list and refresh the per-level uLod uniform */
    prepare: (lodCell: ValueCell<unknown> | undefined, capacity: number) => LodLevelsValue | undefined
    /** cull against the instance grid and append draw entries per lod level */
    cull: (values: CullValues, lodLevels: LodLevelsValue | undefined, cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null, stats: WebGLStats) => void
    /** append the draw entry for the lod level matching the given distance */
    cullSimple: (values: CullValues, lodLevels: LodLevelsValue, d: number, radius: number, scale: number) => void
}

function createMdbList(): MdbList {
    const list: MultiDrawBaseData[] = [];
    let lodLevelsVersion = -1;

    const s = Sphere3D();

    // flattened lod tuples, avoids nested array loads in the hot cell loop
    const lodMin: number[] = [];
    const lodMax: number[] = [];
    const lodCount: number[] = [];

    function prepare(lodCell: ValueCell<unknown> | undefined, capacity: number): LodLevelsValue | undefined {
        const lodLevels = lodCell?.ref.value as LodLevelsValue | undefined;
        const hasLodLevels = !!lodLevels && lodLevels.length > 0;
        const levelCount = hasLodLevels ? lodLevels.length : 1;

        list.length = levelCount;
        const uniformsNeedUpdate = hasLodLevels && lodCell!.ref.version !== lodLevelsVersion;
        for (let j = 0; j < levelCount; ++j) {
            list[j] = getMdbData(capacity, list[j]);
            list[j].count = 0;
            if (hasLodLevels) {
                const l = list[j];
                let created = false;
                if (l.uniforms.length !== 1) {
                    l.uniforms.length = 1;
                    l.uniforms[0] = ['uLod', ValueCell.create(Vec4())];
                    created = true;
                }
                if (uniformsNeedUpdate || created) {
                    ValueCell.update(l.uniforms[0][1], Vec4.set(l.uniforms[0][1].ref.value as Vec4, lodLevels[j][0], lodLevels[j][1], lodLevels[j][2], lodLevels[j][4]));
                }
            } else {
                list[j].uniforms.length = 0;
            }
        }
        if (uniformsNeedUpdate) lodLevelsVersion = lodCell!.ref.version;

        return hasLodLevels ? lodLevels : undefined;
    }

    function emit(l: MultiDrawBaseData, count: number, instanceCount: number, baseInstance: number) {
        const o = l.count;
        if (o > 0 && l.firsts[o - 1] === 0 && l.offsets[o - 1] === 0 && l.counts[o - 1] === count && l.baseInstances[o - 1] + l.instanceCounts[o - 1] === baseInstance) {
            l.instanceCounts[o - 1] += instanceCount;
        } else {
            l.firsts[o] = 0;
            l.offsets[o] = 0;
            l.counts[o] = count;
            l.instanceCounts[o] = instanceCount;
            l.baseInstances[o] = baseInstance;
            l.count += 1;
        }
    }

    function cull(values: CullValues, lodLevels: LodLevelsValue | undefined, cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null, stats: WebGLStats) {
        const drawCount = values.drawCount.ref.value;
        const instanceCount = values.instanceCount.ref.value;
        if (drawCount === 0 || instanceCount === 0) return;

        const [minDistance, maxDistance] = values.uLod.ref.value;
        const hasLod = minDistance !== 0 || maxDistance !== 0;
        const grid = values.instanceGrid.ref.value;

        const levelCount = lodLevels ? lodLevels.length : 0;
        for (let j = 0; j < levelCount; ++j) {
            lodMin[j] = lodLevels![j][0];
            lodMax[j] = lodLevels![j][1];
            lodCount[j] = lodLevels![j][3];
        }

        // caller only invokes this when the instance grid is usable (cellSize > 1)
        const { cellOffsets, cellSpheres, batchOffsets, batchSpheres, batchCount, batchSize } = grid;
        const checkCellOccludedDistance = 2 * batchSize;

        for (let k = 0; k < batchCount; ++k) {
            // batchCell is serial after the grid reorder, cells of a batch are consecutive
            const cBegin = batchOffsets[k];
            const cEnd = batchOffsets[k + 1];
            const cCount = cEnd - cBegin;
            if (cCount === 0) continue;

            s3fromArray(s, batchSpheres, k * 4);
            const bd = p3distanceToPoint(cameraPlane, s.center);
            if (hasLod && (bd + s.radius < minDistance || bd - s.radius > maxDistance)) {
                if (isTimingMode) {
                    stats.culled.lod += cellOffsets[cEnd] - cellOffsets[cBegin];
                }
                continue;
            }
            if (!f3intersectsSphere3D(frustum, s)) {
                if (isTimingMode) {
                    stats.culled.frustum += cellOffsets[cEnd] - cellOffsets[cBegin];
                }
                continue;
            }
            if (isOccluded !== null && isOccluded(s)) {
                if (isTimingMode) {
                    stats.culled.occlusion += cellOffsets[cEnd] - cellOffsets[cBegin];
                }
                continue;
            }

            if (cCount === 1) {
                // single-cell batch: cell sphere equals batch sphere by construction
                // (see calcTopGrid), re-testing it per-cell would be redundant
                const begin = cellOffsets[cBegin];
                const count = cellOffsets[cBegin + 1] - begin;
                if (count === 0) continue;

                if (lodLevels) {
                    const dMax = bd + s.radius;
                    const dMin = bd - s.radius;
                    for (let j = 0; j < levelCount; ++j) {
                        if (dMax < lodMin[j] || dMin > lodMax[j]) continue;
                        emit(list[j], lodCount[j], count, begin);
                    }
                } else {
                    emit(list[0], drawCount, count, begin);
                }
                continue;
            }

            // cells of a fully contained batch don't need frustum tests
            const batchContained = f3containsSphere3D(frustum, s);

            for (let c = cBegin; c < cEnd; ++c) {
                const begin = cellOffsets[c];
                const end = cellOffsets[c + 1];
                const count = end - begin;
                if (count === 0) continue;

                s3fromArray(s, cellSpheres, c * 4);
                const d = p3distanceToPoint(cameraPlane, s.center);
                if (hasLod && (d + s.radius < minDistance || d - s.radius > maxDistance)) {
                    if (isTimingMode) stats.culled.lod += count;
                    continue;
                }
                if (!batchContained && !f3intersectsSphere3D(frustum, s)) {
                    if (isTimingMode) stats.culled.frustum += count;
                    continue;
                }
                if (isOccluded !== null && d - s.radius < checkCellOccludedDistance && isOccluded(s)) {
                    if (isTimingMode) stats.culled.occlusion += count;
                    continue;
                }

                if (lodLevels) {
                    const dMax = d + s.radius;
                    const dMin = d - s.radius;
                    for (let j = 0; j < levelCount; ++j) {
                        if (dMax < lodMin[j] || dMin > lodMax[j]) continue;
                        emit(list[j], lodCount[j], count, begin);
                    }
                } else {
                    emit(list[0], drawCount, count, begin);
                }
            }
        }
    }

    function cullSimple(values: CullValues, lodLevels: LodLevelsValue, d: number, radius: number, scale: number) {
        if (values.drawCount.ref.value === 0 || values.instanceCount.ref.value === 0) return;

        for (let j = 0, jl = lodLevels.length; j < jl; ++j) {
            if (d + radius < lodLevels[j][1] * scale) {
                emit(list[j], lodLevels[j][3], values.instanceCount.ref.value, 0);
                break;
            }
        }
    }

    return { list, prepare, cull, cullSimple };
}

type GraphicsRenderableValues = RenderableValues & BaseValues

export function createRenderable<T extends GraphicsRenderableValues>(renderItem: GraphicsRenderItem, values: T, state: RenderableState): Renderable<T> {
    const id = getNextRenderableId();

    const mdb = createMdbList();
    let cullEnabled = false;
    let lastCullFrame: Frame | undefined;

    return {
        id,
        materialId: renderItem.materialId,
        values,
        state,

        cull: (cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null, stats: WebGLStats, frame: Frame) => {
            // skip recomputation if this renderable was already culled for the current frame
            if (frame === lastCullFrame) {
                stats.cacheHits.cull++;
                return;
            }

            cullEnabled = false;

            if (values.drawCount.ref.value === 0) return;
            if (values.instanceCount.ref.value === 0) return;
            if (values.instanceGrid.ref.value.cellSize <= 1) return;

            const lodLevels = mdb.prepare(values.lodLevels, values.instanceGrid.ref.value.cellCount);
            mdb.cull(values, lodLevels, cameraPlane, frustum, isOccluded, stats);

            cullEnabled = true;
            lastCullFrame = frame;
        },
        uncull: () => {
            cullEnabled = false;
            lastCullFrame = undefined;
        },
        cullSimple: (d: number, radius: number, scale: number) => {
            cullEnabled = false;

            const lodLevels = mdb.prepare(values.lodLevels, Math.max(1, values.instanceGrid.ref.value.cellCount));
            if (!lodLevels) return;

            mdb.cullSimple(values, lodLevels, d, radius, scale);

            cullEnabled = true;
        },
        render: (variant: GraphicsRenderVariant, sharedTexturesCount: number) => {
            if (values.uAlpha && values.alpha) {
                ValueCell.updateIfChanged(values.uAlpha, clamp(values.alpha.ref.value * state.alphaFactor, 0, 1));
            }
            renderItem.render(variant, sharedTexturesCount, cullEnabled ? mdb.list : undefined);
        },
        getByteCount: () => renderItem.getByteCount(),
        getProgram: (variant: GraphicsRenderVariant) => renderItem.getProgram(variant),
        setTransparency: (transparency: Transparency) => renderItem.setTransparency(transparency),
        update: () => {
            renderItem.update();
            lastCullFrame = undefined;
        },
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

        render: () => {
            renderItem.getProgram('compute').finalize(true);
            renderItem.render('compute', 0);
        },
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    };
}