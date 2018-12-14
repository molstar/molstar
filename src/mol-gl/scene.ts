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
        const r = renderables[i]
        if (!r.boundingSphere.radius) continue;
        boundaryHelper.boundaryStep(r.boundingSphere.center, r.boundingSphere.radius);
    }
    boundaryHelper.finishBoundaryStep();
    for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i]
        if (!r.boundingSphere.radius) continue;
        boundaryHelper.extendStep(r.boundingSphere.center, r.boundingSphere.radius);
    }

    Vec3.copy(boundingSphere.center, boundaryHelper.center);
    boundingSphere.radius = boundaryHelper.radius;

    return boundingSphere;
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
    forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
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