/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { Viewport } from 'mol-view/camera/util';
import { Camera } from 'mol-view/camera/base';

import * as glContext from './context'
import Scene, { RenderObject } from './scene';

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
    clear: () => void
    draw: () => void

    setViewport: (viewport: Viewport) => void

    stats: RendererStats
    dispose: () => void
}

const extensions = [
    'OES_element_index_uint',
    'ANGLE_instanced_arrays'
]
const optionalExtensions = [
    'EXT_disjoint_timer_query'
]

namespace Renderer {
    export function create(canvas: HTMLCanvasElement, camera: Camera): Renderer {
        const regl = glContext.create({ canvas, extensions, optionalExtensions, profile: false })
        const scene = Scene.create(regl)

        const baseContext = regl({
            context: {
                model: Mat4.identity(),
                transform: Mat4.identity(),
                view: camera.view,
                projection: camera.projection,
            },
            uniforms: {
                model: regl.context('model' as any),
                transform: regl.context('transform' as any),
                view: regl.context('view' as any),
                projection: regl.context('projection' as any),
                'light.position': Vec3.create(0, 0, -100),
                'light.color': Vec3.create(1.0, 1.0, 1.0),
                'light.ambient': Vec3.create(0.5, 0.5, 0.5),
                'light.falloff': 0,
                'light.radius': 500
            }
        })

        const draw = () => {
            regl.poll() // updates timers and viewport
            camera.update()
            baseContext(state => {
                regl.clear({ color: [0, 0, 0, 1] })
                // TODO painters sort, filter visible, filter picking, visibility culling?
                scene.forEach(r => {
                    r.draw()
                })
            })
        }

        // TODO animate, draw, requestDraw
        return {
            add: (o: RenderObject) => {
                scene.add(o)
            },
            remove: (o: RenderObject) => {
                scene.remove(o)
            },
            clear: () => {
                scene.clear()
            },
            draw,
            setViewport: (viewport: Viewport) => {
                regl({ viewport })
            },
            get stats() {
                return {
                    elementsCount: regl.stats.elementsCount,
                    bufferCount: regl.stats.bufferCount,
                    textureCount: regl.stats.textureCount,
                    shaderCount: regl.stats.shaderCount,
                    renderableCount: scene.count
                }
            },
            dispose: () => {
                regl.destroy()
            }
        }
    }
}

export default Renderer