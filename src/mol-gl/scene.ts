/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { WebGLContext } from './webgl/context';
import { GraphicsRenderObject, createRenderable } from './render-object';
import { Object3D } from './object3d';
import { Sphere3D } from '../mol-math/geometry';
import { CommitQueue } from './commit-queue';
import { now } from '../mol-util/now';
import { arraySetRemove } from '../mol-util/array';
import { BoundaryHelper } from '../mol-math/geometry/boundary-helper';
import { hash1 } from '../mol-data/util';
import { GraphicsRenderable } from './renderable';
import { GraphicsRenderVariants } from './webgl/render-item';
import { clamp } from '../mol-math/interpolate';

const boundaryHelper = new BoundaryHelper('98');

function calculateBoundingSphere(renderables: GraphicsRenderable[], boundingSphere: Sphere3D, onlyVisible: boolean): Sphere3D {
    boundaryHelper.reset();

    for (let i = 0, il = renderables.length; i < il; ++i) {
        if (onlyVisible && !renderables[i].state.visible) continue;

        const boundingSphere = renderables[i].values.boundingSphere.ref.value;
        if (!boundingSphere.radius) continue;

        boundaryHelper.includeSphere(boundingSphere);
    }
    boundaryHelper.finishedIncludeStep();
    for (let i = 0, il = renderables.length; i < il; ++i) {
        if (onlyVisible && !renderables[i].state.visible) continue;

        const boundingSphere = renderables[i].values.boundingSphere.ref.value;
        if (!boundingSphere.radius) continue;

        boundaryHelper.radiusSphere(boundingSphere);
    }

    return boundaryHelper.getSphere(boundingSphere);
}

function renderableSort(a: GraphicsRenderable, b: GraphicsRenderable) {
    const drawProgramIdA = (a.getProgram('colorBlended') || a.getProgram('colorWboit') || a.getProgram('colorDpoit')).id;
    const drawProgramIdB = (b.getProgram('colorBlended') || b.getProgram('colorWboit') || b.getProgram('colorDpoit')).id;
    const materialIdA = a.materialId;
    const materialIdB = b.materialId;

    if (drawProgramIdA !== drawProgramIdB) {
        // sort by program id to minimize gl state changes
        return drawProgramIdA - drawProgramIdB;
    } else if (materialIdA !== materialIdB) {
        // sort by material id to minimize gl state changes
        return materialIdA - materialIdB;
    } else {
        return a.id - b.id;
    }
}

interface Scene extends Object3D {
    readonly count: number
    readonly renderables: ReadonlyArray<GraphicsRenderable>
    readonly boundingSphere: Sphere3D
    readonly boundingSphereVisible: Sphere3D

    readonly primitives: Scene.Group
    readonly volumes: Scene.Group

    /** Returns `true` if some visibility has changed, `false` otherwise. */
    syncVisibility: () => boolean
    update: (objects: ArrayLike<GraphicsRenderObject> | undefined, keepBoundingSphere: boolean) => void
    add: (o: GraphicsRenderObject) => void // GraphicsRenderable
    remove: (o: GraphicsRenderObject) => void
    commit: (maxTimeMs?: number) => boolean
    readonly needsCommit: boolean
    has: (o: GraphicsRenderObject) => boolean
    clear: () => void
    forEach: (callbackFn: (value: GraphicsRenderable, key: GraphicsRenderObject) => void) => void
    /** Marker average of primitive renderables */
    readonly markerAverage: number
    /** Opacity average of primitive renderables */
    readonly opacityAverage: number
    /** Is `true` if any primitive renderable (possibly) has any opaque part */
    readonly hasOpaque: boolean
}

namespace Scene {
    export interface Group extends Object3D {
        readonly renderables: ReadonlyArray<GraphicsRenderable>
    }

    export function create(ctx: WebGLContext, variants = GraphicsRenderVariants): Scene {
        const renderableMap = new Map<GraphicsRenderObject, GraphicsRenderable>();
        const renderables: GraphicsRenderable[] = [];
        const boundingSphere = Sphere3D();
        const boundingSphereVisible = Sphere3D();

        const primitives: GraphicsRenderable[] = [];
        const volumes: GraphicsRenderable[] = [];

        let boundingSphereDirty = true;
        let boundingSphereVisibleDirty = true;

        let markerAverage = 0;
        let opacityAverage = 0;
        let hasOpaque = false;

        const object3d = Object3D.create();
        const { view, position, direction, up } = object3d;

        function add(o: GraphicsRenderObject) {
            if (!renderableMap.has(o)) {
                const renderable = createRenderable(ctx, o, variants);
                renderables.push(renderable);
                if (o.type === 'direct-volume') {
                    volumes.push(renderable);
                } else {
                    primitives.push(renderable);
                }
                renderableMap.set(o, renderable);
                boundingSphereDirty = true;
                boundingSphereVisibleDirty = true;
                return renderable;
            } else {
                console.warn(`RenderObject with id '${o.id}' already present`);
                return renderableMap.get(o)!;
            }
        }

        function remove(o: GraphicsRenderObject) {
            const renderable = renderableMap.get(o);
            if (renderable) {
                renderable.dispose();
                arraySetRemove(renderables, renderable);
                arraySetRemove(primitives, renderable);
                arraySetRemove(volumes, renderable);
                renderableMap.delete(o);
                boundingSphereDirty = true;
                boundingSphereVisibleDirty = true;
            }
        }

        const commitBulkSize = 100;
        function commit(maxTimeMs: number) {
            const start = now();

            let i = 0;

            while (true) {
                const o = commitQueue.tryGetRemove();
                if (!o) break;
                remove(o);
                if (++i % commitBulkSize === 0 && now() - start > maxTimeMs) return false;
            }

            while (true) {
                const o = commitQueue.tryGetAdd();
                if (!o) break;
                add(o);
                if (++i % commitBulkSize === 0 && now() - start > maxTimeMs) return false;
            }

            renderables.sort(renderableSort);
            markerAverage = calculateMarkerAverage();
            opacityAverage = calculateOpacityAverage();
            hasOpaque = calculateHasOpaque();
            return true;
        }

        const commitQueue = new CommitQueue();

        let visibleHash = -1;
        function computeVisibleHash() {
            let hash = 23;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                if (!renderables[i].state.visible) continue;
                hash = (31 * hash + renderables[i].id) | 0;
            }
            hash = hash1(hash);
            if (hash === -1) hash = 0;
            return hash;
        }

        function syncVisibility() {
            const newVisibleHash = computeVisibleHash();
            if (newVisibleHash !== visibleHash) {
                boundingSphereVisibleDirty = true;
                markerAverage = calculateMarkerAverage();
                opacityAverage = calculateOpacityAverage();
                hasOpaque = calculateHasOpaque();
                visibleHash = newVisibleHash;
                return true;
            } else {
                return false;
            }
        }

        function calculateMarkerAverage() {
            if (primitives.length === 0) return 0;
            let count = 0;
            let markerAverage = 0;
            for (let i = 0, il = primitives.length; i < il; ++i) {
                if (!primitives[i].state.visible) continue;
                markerAverage += primitives[i].values.markerAverage.ref.value;
                count += 1;
            }
            return count > 0 ? markerAverage / count : 0;
        }

        function calculateOpacityAverage() {
            if (primitives.length === 0) return 0;
            let count = 0;
            let opacityAverage = 0;
            for (let i = 0, il = primitives.length; i < il; ++i) {
                const p = primitives[i];
                if (!p.state.visible) continue;
                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(p.values.alpha.ref.value * p.state.alphaFactor, 0, 1);
                const xray = p.values.dXrayShaded?.ref.value ? 0.5 : 1;
                const fuzzy = p.values.dPointStyle?.ref.value === 'fuzzy' ? 0.5 : 1;
                const text = p.values.dGeometryType.ref.value === 'text' ? 0.5 : 1;
                opacityAverage += (1 - p.values.transparencyAverage.ref.value) * alpha * xray * fuzzy * text;
                count += 1;
            }
            return count > 0 ? opacityAverage / count : 0;
        }

        function calculateHasOpaque() {
            if (primitives.length === 0) return false;
            for (let i = 0, il = primitives.length; i < il; ++i) {
                const p = primitives[i];
                if (!p.state.visible) continue;

                if (p.state.opaque) return true;
                if (p.state.alphaFactor === 1 && p.values.alpha.ref.value === 1 && p.values.transparencyAverage.ref.value !== 1) return true;
                if (p.values.dTransparentBackfaces?.ref.value === 'opaque') return true;
            }
            return false;
        }

        return {
            view, position, direction, up,

            renderables,
            primitives: { view, position, direction, up, renderables: primitives },
            volumes: { view, position, direction, up, renderables: volumes },

            syncVisibility,
            update(objects, keepBoundingSphere) {
                Object3D.update(object3d);
                if (objects) {
                    for (let i = 0, il = objects.length; i < il; ++i) {
                        renderableMap.get(objects[i])?.update();
                    }
                } else {
                    for (let i = 0, il = renderables.length; i < il; ++i) {
                        renderables[i].update();
                    }
                }
                if (!keepBoundingSphere) {
                    boundingSphereDirty = true;
                    boundingSphereVisibleDirty = true;
                } else {
                    syncVisibility();
                }
                markerAverage = calculateMarkerAverage();
                opacityAverage = calculateOpacityAverage();
                hasOpaque = calculateHasOpaque();
            },
            add: (o: GraphicsRenderObject) => commitQueue.add(o),
            remove: (o: GraphicsRenderObject) => commitQueue.remove(o),
            commit: (maxTime = Number.MAX_VALUE) => commit(maxTime),
            get needsCommit() { return !commitQueue.isEmpty; },
            has: (o: GraphicsRenderObject) => {
                return renderableMap.has(o);
            },
            clear: () => {
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    renderables[i].dispose();
                }
                renderables.length = 0;
                primitives.length = 0;
                volumes.length = 0;
                renderableMap.clear();
                boundingSphereDirty = true;
                boundingSphereVisibleDirty = true;
            },
            forEach: (callbackFn: (value: GraphicsRenderable, key: GraphicsRenderObject) => void) => {
                renderableMap.forEach(callbackFn);
            },
            get count() {
                return renderables.length;
            },
            get boundingSphere() {
                if (boundingSphereDirty) {
                    calculateBoundingSphere(renderables, boundingSphere, false);
                    boundingSphereDirty = false;
                }
                return boundingSphere;
            },
            get boundingSphereVisible() {
                if (boundingSphereVisibleDirty) {
                    calculateBoundingSphere(renderables, boundingSphereVisible, true);
                    boundingSphereVisibleDirty = false;
                }
                return boundingSphereVisible;
            },
            get markerAverage() {
                return markerAverage;
            },
            get opacityAverage() {
                return opacityAverage;
            },
            get hasOpaque() {
                return hasOpaque;
            },
        };
    }
}

export { Scene };