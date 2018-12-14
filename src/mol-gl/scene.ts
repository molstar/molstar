/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable } from './renderable'
import { WebGLContext } from './webgl/context';
import { RenderableValues, BaseValues } from './renderable/schema';
import { RenderObject, createRenderable } from './render-object';
import { Object3D } from './object3d';
import { Sphere3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import { BoundaryHelper } from 'mol-math/geometry/boundary-helper';

const boundaryHelper = new BoundaryHelper();
function calculateBoundingSphere(renderables: Renderable<RenderableValues & BaseValues>[], boundingSphere: Sphere3D): Sphere3D {
    boundaryHelper.reset(0.1);

    for (let i = 0, il = renderables.length; i < il; ++i) {
        const boundingSphere = renderables[i].values.boundingSphere.ref.value
        if (!boundingSphere.radius) continue;
        boundaryHelper.boundaryStep(boundingSphere.center, boundingSphere.radius);
    }
    boundaryHelper.finishBoundaryStep();
    for (let i = 0, il = renderables.length; i < il; ++i) {
        const boundingSphere = renderables[i].values.boundingSphere.ref.value
        if (!boundingSphere.radius) continue;
        boundaryHelper.extendStep(boundingSphere.center, boundingSphere.radius);
    }

    Vec3.copy(boundingSphere.center, boundaryHelper.center);
    boundingSphere.radius = boundaryHelper.radius;

    return boundingSphere;
}

function renderableSort(a: Renderable<RenderableValues & BaseValues>, b: Renderable<RenderableValues & BaseValues>) {
    const drawProgramIdA = a.getProgram('draw').id
    const drawProgramIdB = b.getProgram('draw').id
    const zA = a.values.boundingSphere.ref.value.center[2]
    const zB = a.values.boundingSphere.ref.value.center[2]

    if (drawProgramIdA !== drawProgramIdB) {
        return drawProgramIdA - drawProgramIdB; // sort by program id to minimize gl state changes
    } else if (zA !== zB) {
        return a.state.opaque
            ? zA - zB // when opaque draw closer elements first to minimize overdraw
            : zB - zA // when transparent draw elements last to maximize partial visibility
    } else {
        return a.id - b.id;
    }
}

interface Scene extends Object3D {
    readonly count: number
    readonly renderables: ReadonlyArray<Renderable<RenderableValues & BaseValues>>
    readonly boundingSphere: Sphere3D

    update: (keepBoundingSphere?: boolean) => void
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    has: (o: RenderObject) => boolean
    clear: () => void
    forEach: (callbackFn: (value: Renderable<RenderableValues & BaseValues>, key: RenderObject) => void) => void
}

namespace Scene {
    export function create(ctx: WebGLContext): Scene {
        const renderableMap = new Map<RenderObject, Renderable<RenderableValues & BaseValues>>()
        const renderables: Renderable<RenderableValues & BaseValues>[] = []
        const boundingSphere = Sphere3D.zero()
        let boundingSphereDirty = true

        const object3d = Object3D.create()

        return {
            get view () { return object3d.view },
            get position () { return object3d.position },
            get direction () { return object3d.direction },
            get up () { return object3d.up },

            update: (keepBoundingSphere?: boolean) => {
                Object3D.update(object3d)
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    renderables[i].update()
                }
                if (!keepBoundingSphere) boundingSphereDirty = true
            },
            add: (o: RenderObject) => {
                if (!renderableMap.has(o)) {
                    const renderable = createRenderable(ctx, o)
                    renderables.push(renderable)
                    renderables.sort(renderableSort)
                    renderableMap.set(o, renderable)
                    boundingSphereDirty = true
                } else {
                    console.warn(`RenderObject with id '${o.id}' already present`)
                }
            },
            remove: (o: RenderObject) => {
                const renderable = renderableMap.get(o)
                if (renderable) {
                    renderable.dispose()
                    renderables.splice(renderables.indexOf(renderable), 1)
                    renderables.sort(renderableSort)
                    renderableMap.delete(o)
                    boundingSphereDirty = true
                }
            },
            has: (o: RenderObject) => {
                return renderableMap.has(o)
            },
            clear: () => {
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    renderables[i].dispose()
                }
                renderables.length = 0
                renderableMap.clear()
                boundingSphereDirty = true
            },
            forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach(callbackFn)
            },
            get count() {
                return renderables.length
            },
            renderables,
            get boundingSphere() {
                if (boundingSphereDirty) calculateBoundingSphere(renderables, boundingSphere)
                boundingSphereDirty = false
                return boundingSphere
            }
        }
    }
}

export default Scene