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

function calculateBoundingSphere(renderableMap: Map<RenderObject, Renderable<RenderableValues & BaseValues>>, boundingSphere: Sphere3D): Sphere3D {
    // let count = 0
    // const center = Vec3.set(boundingSphere.center, 0, 0, 0)
    // renderableMap.forEach(r => {
    //     if (r.boundingSphere.radius) {
    //         Vec3.add(center, center, r.boundingSphere.center)
    //         ++count
    //     }
    // })
    // if (count > 0) {
    //     Vec3.scale(center, center, 1 / count)
    // }

    // let radius = 0
    // renderableMap.forEach(r => {
    //     if (r.boundingSphere.radius) {
    //         radius = Math.max(radius, Vec3.distance(center, r.boundingSphere.center) + r.boundingSphere.radius)
    //     }
    // })
    // boundingSphere.radius = radius

    const spheres: Sphere3D[] = [];
    renderableMap.forEach(r => {
        if (!r.state.visible || !r.boundingSphere.radius) return;
        spheres.push(r.boundingSphere)
    });
    const bs = Sphere3D.getBoundingSphereFromSpheres(spheres, 0.1);

    Vec3.copy(boundingSphere.center, bs.center);
    boundingSphere.radius = bs.radius;

    return boundingSphere;
}

interface Scene extends Object3D {
    readonly count: number
    readonly boundingSphere: Sphere3D

    update: () => void
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    clear: () => void
    forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
    eachOpaque: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
    eachTransparent: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
}

namespace Scene {
    export function create(ctx: WebGLContext): Scene {
        const renderableMap = new Map<RenderObject, Renderable<RenderableValues & BaseValues>>()
        const boundingSphere = Sphere3D.zero()
        let boundingSphereDirty = true

        const object3d = Object3D.create()

        return {
            get view () { return object3d.view },
            get position () { return object3d.position },
            get direction () { return object3d.direction },
            get up () { return object3d.up },

            update: () => {
                Object3D.update(object3d)
                renderableMap.forEach(r => r.update())
                boundingSphereDirty = true
            },

            add: (o: RenderObject) => {
                if (!renderableMap.has(o)) {
                    renderableMap.set(o, createRenderable(ctx, o))
                    boundingSphereDirty = true
                } else {
                    console.warn(`RenderObject with id '${o.id}' already present`)
                }
            },
            remove: (o: RenderObject) => {
                const renderable = renderableMap.get(o)
                if (renderable) {
                    renderable.dispose()
                    renderableMap.delete(o)
                    boundingSphereDirty = true
                }
            },
            clear: () => {
                renderableMap.forEach(renderable => renderable.dispose())
                renderableMap.clear()
                boundingSphereDirty = true
            },
            forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach(callbackFn)
            },
            eachOpaque: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach((r, o) => {
                    if (r.opaque) callbackFn(r, o)
                })
            },
            eachTransparent: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach((r, o) => {
                    if (!r.opaque) callbackFn(r, o)
                })
            },
            get count() {
                return renderableMap.size
            },
            get boundingSphere() {
                if (boundingSphereDirty) calculateBoundingSphere(renderableMap, boundingSphere)
                boundingSphereDirty = false
                return boundingSphere
            }
        }
    }
}

export default Scene