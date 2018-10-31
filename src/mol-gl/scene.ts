/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable } from './renderable'
import { WebGLContext } from './webgl/context';
import { RenderableValues, BaseValues } from './renderable/schema';
import { RenderObject, createRenderable } from './render-object';
import { Object3D, createObject3D } from './object3d';
import { Sphere3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';

function calculateBoundingSphere(renderableMap: Map<RenderObject, Renderable<RenderableValues & BaseValues>>): Sphere3D {
    let count = 0
    const center = Vec3.zero()
    renderableMap.forEach(r => {
        if (r.boundingSphere.radius) {
            Vec3.add(center, center, r.boundingSphere.center)
            ++count
        }
    })
    if (count > 0) {
        Vec3.scale(center, center, 1 / count)
    }

    let radius = 0
    renderableMap.forEach(r => {
        if (r.boundingSphere.radius) {
            radius = Math.max(radius, Vec3.distance(center, r.boundingSphere.center) + r.boundingSphere.radius)
        }
    })

    return { center, radius };
}

interface Scene extends Object3D {
    readonly count: number
    readonly boundingSphere: Sphere3D

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
        let boundingSphere: Sphere3D | undefined

        const { view, position, up, direction, update } = createObject3D()

        return {
            // ...createObject3D(), // TODO does not work in conjunction with getter
            view, position, up, direction,

            update: () => {
                update()
                renderableMap.forEach(r => r.update())
                boundingSphere = undefined
            },

            add: (o: RenderObject) => {
                if (!renderableMap.has(o)) {
                    renderableMap.set(o, createRenderable(ctx, o))
                    boundingSphere = undefined
                } else {
                    console.warn(`RenderObject with id '${o.id}' already present`)
                }
            },
            remove: (o: RenderObject) => {
                const renderable = renderableMap.get(o)
                if (renderable) {
                    renderable.dispose()
                    renderableMap.delete(o)
                    boundingSphere = undefined
                }
            },
            clear: () => {
                renderableMap.forEach(renderable => renderable.dispose())
                renderableMap.clear()
                boundingSphere = undefined
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
                if (boundingSphere) return boundingSphere
                // TODO avoid array creation
                boundingSphere = calculateBoundingSphere(renderableMap)
                return boundingSphere
            }
        }
    }
}

export default Scene