/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { Viewport } from 'mol-view/camera/util';
import { Camera } from 'mol-view/camera/base';

import Scene, { RenderObject } from './scene';
import { Context } from './webgl/context';

export interface RendererStats {
    elementsCount: number
    bufferCount: number
    textureCount: number
    shaderCount: number
    renderableCount: number
}

interface Renderer {
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    update: () => void
    clear: () => void
    draw: () => void

    setViewport: (viewport: Viewport) => void

    stats: RendererStats
    dispose: () => void
}

// function getPixelRatio() {
//     return (typeof window !== 'undefined') ? window.devicePixelRatio : 1
// }

namespace Renderer {
    export function create(ctx: Context, camera: Camera): Renderer {
        const { gl } = ctx
        const scene = Scene.create(ctx)

        // const baseContext = regl({
        //     context: {
        //         model: Mat4.identity(),
        //         transform: Mat4.identity(),
        //         view: camera.view,
        //         projection: camera.projection,
        //     },
        //     uniforms: {
        //         pixelRatio: getPixelRatio(),
        //         viewportHeight: regl.context('viewportHeight'),

        //         model: regl.context('model' as any),
        //         transform: regl.context('transform' as any),
        //         view: regl.context('view' as any),
        //         projection: regl.context('projection' as any),

        //         'light.position': Vec3.create(0, 0, -100),
        //         'light.color': Vec3.create(1.0, 1.0, 1.0),
        //         'light.ambient': Vec3.create(0.5, 0.5, 0.5),
        //         'light.falloff': 0,
        //         'light.radius': 500
        //     }
        // })

        const draw = () => {
            // TODO clear color
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            // TODO painters sort, filter visible, filter picking, visibility culling?
            scene.forEach((r, o) => {
                if (o.visible) r.draw()
            })
        }

        return {
            add: (o: RenderObject) => {
                scene.add(o)
            },
            remove: (o: RenderObject) => {
                scene.remove(o)
            },
            update: () => {
                scene.forEach((r, o) => r.update(o))
            },
            clear: () => {
                scene.clear()
            },
            draw,
            setViewport: (viewport: Viewport) => {
                gl.viewport(viewport.x, viewport.y, viewport.width, viewport.height)
            },
            get stats() {
                return {
                    // elementsCount: regl.stats.elementsCount,
                    // bufferCount: regl.stats.bufferCount,
                    // textureCount: regl.stats.textureCount,
                    // shaderCount: regl.stats.shaderCount,
                    // renderableCount: scene.count
                } as any
            },
            dispose: () => {
                scene.clear()
            }
        }
    }
}

export default Renderer