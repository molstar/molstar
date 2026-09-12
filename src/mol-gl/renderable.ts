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

    cull: (cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null, stats: WebGLStats) => void
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

/**
 * Tracks the inputs used for the last successful cull so a renderable can
 * decide for itself whether anything relevant actually changed, instead of
 * trusting an externally managed "same frame" signal (which is easy to get
 * wrong, e.g. across stereo eyes or multi-sample jitter passes).
 */
interface CullCache {
    /** true if `values`, `cameraPlane`, `frustum` and `isOccluded` are all unchanged since the last `commit` */
    unchanged: (values: CullValues, cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null) => boolean
    /** snapshot the inputs used for a just-completed cull */
    commit: (values: CullValues, cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null) => void
    /** force the next `unchanged` call to return false, e.g. after `values` were updated */
    invalidate: () => void
}

function createCullCache(): CullCache {
    let has = false;
    let drawCountVersion = -1;
    let instanceCountVersion = -1;
    let instanceGridVersion = -1;
    let uLodVersion = -1;
    let lodLevelsVersion = -1;
    let isOccludedRef: ((s: Sphere3D) => boolean) | null = null;
    // cameraPlane: normal.xyz + constant; frustum: 6 planes of the same
    // Vec3/constant are plain (double-precision) numbers, Float32 would round and never compare equal
    const plane = new Float64Array(4);
    const frustumArray = new Float64Array(24);

    // AA jitter perturbs the projection by a sub-pixel offset every sample, so exact
    // equality would never hold and cull would rerun for every jitter sample; a small
    // relative tolerance absorbs that. Measured jitter-only relative delta on `constant`
    // was ~2.2e-3 (JitterVectors offsets go up to ~0.5px, more of a shift for smaller
    // viewports) - 1e-3 was too tight and still triggered a recull every sample; 1e-2
    // gives real margin while still far below any perceptible real camera move.
    const CAMERA_REL_EPSILON = 1e-2;
    function nearlyEqual(a: number, b: number) {
        return Math.abs(a - b) <= CAMERA_REL_EPSILON * Math.max(1, Math.abs(a), Math.abs(b));
    }

    function ownValuesUnchanged(values: CullValues) {
        return values.drawCount.ref.version === drawCountVersion &&
            values.instanceCount.ref.version === instanceCountVersion &&
            values.instanceGrid.ref.version === instanceGridVersion &&
            values.uLod.ref.version === uLodVersion &&
            (values.lodLevels ? values.lodLevels.ref.version : -1) === lodLevelsVersion;
    }

    function cameraUnchanged(cameraPlane: Plane3D, frustum: Frustum3D) {
        const n = cameraPlane.normal;
        if (!nearlyEqual(plane[0], n[0]) || !nearlyEqual(plane[1], n[1]) || !nearlyEqual(plane[2], n[2]) || !nearlyEqual(plane[3], cameraPlane.constant)) return false;
        let o = 0;
        const planes = frustum as unknown as Plane3D[];
        for (let i = 0; i < 6; ++i) {
            const p = planes[i];
            const pn = p.normal;
            if (!nearlyEqual(frustumArray[o], pn[0]) || !nearlyEqual(frustumArray[o + 1], pn[1]) || !nearlyEqual(frustumArray[o + 2], pn[2]) || !nearlyEqual(frustumArray[o + 3], p.constant)) return false;
            o += 4;
        }
        return true;
    }

    return {
        unchanged(values, cameraPlane, frustum, isOccluded) {
            if (!has) return false;
            if (isOccluded !== isOccludedRef) return false;
            // check camera first: for large scenes cameraUnchanged is the one most likely
            // to be the culprit if it never becomes true again after the first camera move
            return cameraUnchanged(cameraPlane, frustum) && ownValuesUnchanged(values);
        },
        commit(values, cameraPlane, frustum, isOccluded) {
            has = true;
            drawCountVersion = values.drawCount.ref.version;
            instanceCountVersion = values.instanceCount.ref.version;
            instanceGridVersion = values.instanceGrid.ref.version;
            uLodVersion = values.uLod.ref.version;
            lodLevelsVersion = values.lodLevels ? values.lodLevels.ref.version : -1;
            isOccludedRef = isOccluded;
            const n = cameraPlane.normal;
            plane[0] = n[0]; plane[1] = n[1]; plane[2] = n[2]; plane[3] = cameraPlane.constant;
            let o = 0;
            const planes = frustum as unknown as Plane3D[];
            for (let i = 0; i < 6; ++i) {
                const p = planes[i];
                const pn = p.normal;
                frustumArray[o] = pn[0]; frustumArray[o + 1] = pn[1]; frustumArray[o + 2] = pn[2]; frustumArray[o + 3] = p.constant;
                o += 4;
            }
        },
        invalidate() {
            has = false;
        },
    };
}

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
    const cullCache = createCullCache();

    return {
        id,
        materialId: renderItem.materialId,
        values,
        state,

        cull: (cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null, stats: WebGLStats) => {
            // skip recomputation if nothing relevant to culling has actually changed
            if (cullCache.unchanged(values, cameraPlane, frustum, isOccluded)) {
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
            cullCache.commit(values, cameraPlane, frustum, isOccluded);
        },
        uncull: () => {
            cullEnabled = false;
            cullCache.invalidate();
        },
        cullSimple: (d: number, radius: number, scale: number) => {
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
            cullCache.invalidate();
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