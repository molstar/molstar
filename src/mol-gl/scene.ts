/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable } from './renderable'
import { Context } from './webgl/context';
import { RenderableValues } from './renderable/schema';
import { RenderObject, createRenderable } from './render-object';


interface Scene {
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    clear: () => void
    forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
    eachOpaque: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
    eachTransparent: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
    count: number
}

namespace Scene {
    export function create(ctx: Context): Scene {
        const renderableMap = new Map<RenderObject, Renderable<RenderableValues>>()

        return {
            add: (o: RenderObject) => {
                if (!renderableMap.has(o)) {
                    renderableMap.set(o, createRenderable(ctx, o))
                } else {
                    console.warn(`RenderObject with id '${o.id}' already present`)
                }
            },
            remove: (o: RenderObject) => {
                const renderable = renderableMap.get(o)
                if (renderable) {
                    renderable.dispose()
                    renderableMap.delete(o)
                }
            },
            clear: () => {
                renderableMap.forEach(renderable => renderable.dispose())
                renderableMap.clear()
            },
            forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach(callbackFn)
            },
            eachOpaque: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach((r, o) => {
                    if (o.values.uAlpha.ref.value === 1) callbackFn(r, o)
                })
            },
            eachTransparent: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach((r, o) => {
                    if (o.values.uAlpha.ref.value < 1) callbackFn(r, o)
                })
            },
            get count() {
                return renderableMap.size
            }
        }
    }
}

export default Scene