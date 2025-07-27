/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { WebGLContext } from './webgl/context';
import { GraphicsRenderObject, createRenderable } from './render-object';
import { Object3D } from './object3d';
import { Sphere3D } from '../mol-math/geometry/primitives/sphere3d';
import { CommitQueue } from './commit-queue';
import { now } from '../mol-util/now';
import { arraySetRemove } from '../mol-util/array';
import { BoundaryHelper } from '../mol-math/geometry/boundary-helper';
import { hash1 } from '../mol-data/util';
import { GraphicsRenderable } from './renderable';
import { Transparency } from './webgl/render-item';
import { clamp } from '../mol-math/interpolate';

const boundaryHelper = new BoundaryHelper('98');

function calculateBoundingSphere(renderables: GraphicsRenderable[], boundingSphere: Sphere3D, onlyVisible: boolean): Sphere3D {
    boundaryHelper.reset();

    for (let i = 0, il = renderables.length; i < il; ++i) {
        if (onlyVisible && !renderables[i].state.visible) continue;
        if (!renderables[i].values.drawCount.ref.value) continue;

        const boundingSphere = renderables[i].values.boundingSphere.ref.value;
        if (!boundingSphere.radius) continue;

        boundaryHelper.includeSphere(boundingSphere);
    }
    boundaryHelper.finishedIncludeStep();
    for (let i = 0, il = renderables.length; i < il; ++i) {
        if (onlyVisible && !renderables[i].state.visible) continue;
        if (!renderables[i].values.drawCount.ref.value) continue;

        const boundingSphere = renderables[i].values.boundingSphere.ref.value;
        if (!boundingSphere.radius) continue;

        boundaryHelper.radiusSphere(boundingSphere);
    }

    return boundaryHelper.getSphere(boundingSphere);
}

function renderableSort(a: GraphicsRenderable, b: GraphicsRenderable) {
    const drawProgramIdA = a.getProgram('color').id;
    const drawProgramIdB = b.getProgram('color').id;
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
    setTransparency: (transparency: Transparency) => void
    update: (objects: ArrayLike<GraphicsRenderObject> | undefined, keepBoundingSphere: boolean) => void
    add: (o: GraphicsRenderObject) => void // GraphicsRenderable
    remove: (o: GraphicsRenderObject) => void
    commit: (maxTimeMs?: number) => boolean
    readonly needsCommit: boolean
    readonly commitQueueSize: number
    has: (o: GraphicsRenderObject) => boolean
    clear: () => void
    forEach: (callbackFn: (value: GraphicsRenderable, key: GraphicsRenderObject) => void) => void
    /** Marker average of primitive renderables */
    readonly markerAverage: number
    /** Emissive average of primitive renderables */
    readonly emissiveAverage: number
    /** Opacity average of primitive renderables */
    readonly opacityAverage: number
    /** Transparency minimum, excluding fully opaque, of primitive renderables */
    readonly transparencyMin: number
    /** Is `true` if any primitive renderable (possibly) has any opaque part */
    readonly hasOpaque: boolean
}

namespace Scene {
    export interface Group extends Object3D {
        readonly renderables: ReadonlyArray<GraphicsRenderable>
    }

    export function create(ctx: WebGLContext, transparency: Transparency = 'blended'): Scene {
        const renderableMap = new Map<GraphicsRenderObject, GraphicsRenderable>();
        const renderables: GraphicsRenderable[] = [];
        const boundingSphere = Sphere3D();
        const boundingSphereVisible = Sphere3D();

        const primitives: GraphicsRenderable[] = [];
        const volumes: GraphicsRenderable[] = [];

        let boundingSphereDirty = true;
        let boundingSphereVisibleDirty = true;

        let markerAverageDirty = true;
        let emissiveAverageDirty = true;
        let opacityAverageDirty = true;
        let transparencyMinDirty = true;
        let hasOpaqueDirty = true;

        let markerAverage = 0;
        let emissiveAverage = 0;
        let opacityAverage = 0;
        let transparencyMin = 0;
        let hasOpaque = false;

        const object3d = Object3D.create();
        const { view, position, direction, up } = object3d;

        function add(o: GraphicsRenderObject) {
            if (!renderableMap.has(o)) {
                const renderable = createRenderable(ctx, o, transparency);
                renderables.push(renderable);
                if (o.type === 'direct-volume') {
                    volumes.push(renderable);
                } else {
                    primitives.push(renderable);
                }
                renderableMap.set(o, renderable);
                boundingSphereDirty = true;
                boundingSphereVisibleDirty = true;
            } else {
                console.warn(`RenderObject with id '${o.id}' already present`);
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
            markerAverageDirty = true;
            emissiveAverageDirty = true;
            opacityAverageDirty = true;
            transparencyMinDirty = true;
            hasOpaqueDirty = true;
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
                markerAverageDirty = true;
                emissiveAverageDirty = true;
                opacityAverageDirty = true;
                transparencyMinDirty = true;
                hasOpaqueDirty = true;
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

        function calculateEmissiveAverage() {
            if (primitives.length === 0) return 0;
            let count = 0;
            let emissiveAverage = 0;
            for (let i = 0, il = primitives.length; i < il; ++i) {
                if (!primitives[i].state.visible) continue;
                emissiveAverage += primitives[i].values.emissiveAverage.ref.value + primitives[i].values.uEmissive.ref.value;
                count += 1;
            }
            return count > 0 ? emissiveAverage / count : 0;
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
                const xray = (p.values.dXrayShaded?.ref.value === 'on' || p.values.dXrayShaded?.ref.value === 'inverted') ? 0.5 : 1;
                const fuzzy = p.values.dPointStyle?.ref.value === 'fuzzy' ? 0.5 : 1;
                const text = p.values.dGeometryType.ref.value === 'text' ? 0.5 : 1;
                const image = p.values.dGeometryType.ref.value === 'image' ? 0.5 : 1;
                opacityAverage += (1 - p.values.transparencyAverage.ref.value) * alpha * xray * fuzzy * text * image;
                count += 1;
            }
            return count > 0 ? opacityAverage / count : 0;
        }

        /** exclude fully opaque parts */
        function calculateTransparencyMin() {
            if (primitives.length === 0) return 1;
            let transparencyMin = 1;
            const transparenyValues: number[] = [];
            for (let i = 0, il = primitives.length; i < il; ++i) {
                const p = primitives[i];
                if (!p.state.visible) continue;
                transparenyValues.length = 0;
                const alpha = clamp(p.values.alpha.ref.value * p.state.alphaFactor, 0, 1);
                if (alpha < 1) transparenyValues.push(1 - alpha);
                if (p.values.dXrayShaded?.ref.value === 'on' ||
                    p.values.dXrayShaded?.ref.value === 'inverted' ||
                    p.values.dPointStyle?.ref.value === 'fuzzy' ||
                    p.values.dGeometryType.ref.value === 'text' ||
                    p.values.dGeometryType.ref.value === 'image'
                ) transparenyValues.push(0.5);
                if (p.values.transparencyMin.ref.value > 0) transparenyValues.push(p.values.transparencyMin.ref.value);
                transparencyMin = Math.min(transparencyMin, ...transparenyValues);
            }
            return transparencyMin;
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
            setTransparency: (value: Transparency) => {
                transparency = value;
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    renderables[i].setTransparency(value);
                }
            },
            update(objects, keepBoundingSphere) {
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
                markerAverageDirty = true;
                emissiveAverageDirty = true;
                opacityAverageDirty = true;
                transparencyMinDirty = true;
                hasOpaqueDirty = true;
            },
            add: (o: GraphicsRenderObject) => commitQueue.add(o),
            remove: (o: GraphicsRenderObject) => commitQueue.remove(o),
            commit: (maxTime = Number.MAX_VALUE) => commit(maxTime),
            get commitQueueSize() { return commitQueue.size; },
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
                if (markerAverageDirty) {
                    markerAverage = calculateMarkerAverage();
                    markerAverageDirty = false;
                }
                return markerAverage;
            },
            get emissiveAverage() {
                if (emissiveAverageDirty) {
                    emissiveAverage = calculateEmissiveAverage();
                    emissiveAverageDirty = false;
                }
                return emissiveAverage;
            },
            get opacityAverage() {
                if (opacityAverageDirty) {
                    opacityAverage = calculateOpacityAverage();
                    opacityAverageDirty = false;
                }
                return opacityAverage;
            },
            get transparencyMin() {
                if (transparencyMinDirty) {
                    transparencyMin = calculateTransparencyMin();
                    transparencyMinDirty = false;
                }
                return transparencyMin;
            },
            get hasOpaque() {
                if (hasOpaqueDirty) {
                    hasOpaque = calculateHasOpaque();
                    hasOpaqueDirty = false;
                }
                return hasOpaque;
            },
        };
    }
}

export { Scene };