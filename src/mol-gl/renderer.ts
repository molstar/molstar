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
import { Mat4, Vec3 } from 'mol-math/linear-algebra';

export interface RendererStats {
    renderableCount: number
    programCount: number
    shaderCount: number
    bufferCount: number
    textureCount: number
    vaoCount: number
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

function getPixelRatio() {
    return (typeof window !== 'undefined') ? window.devicePixelRatio : 1
}

namespace Renderer {
    export function create(ctx: Context, camera: Camera): Renderer {
        const { gl } = ctx
        const scene = Scene.create(ctx)

        const model = Mat4.identity()
        const viewport = Viewport.create(0, 0, 0, 0)
        const pixelRatio = getPixelRatio()

        // const light_position = Vec3.create(0, 0, -100)
        const light_color = Vec3.create(1.0, 1.0, 1.0)
        const light_ambient = Vec3.create(0.5, 0.5, 0.5)

        const draw = () => {
            // TODO clear color
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
            gl.enable(gl.DEPTH_TEST)

            // TODO painters sort, filter visible, filter picking, visibility culling?
            let currentProgramId = -1
            scene.forEach((r, o) => {
                if (o.visible) {
                    if (currentProgramId !== r.program.id) {
                        r.program.use()
                        r.program.setUniforms({
                            model,
                            view: camera.view,
                            projection: camera.projection,

                            pixelRatio,
                            viewportHeight: viewport.height,

                            // light_position,
                            light_color,
                            light_ambient,
                        })
                        currentProgramId = r.program.id
                    }
                    r.draw()
                }
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
            setViewport: (newViewport: Viewport) => {
                Viewport.copy(viewport, newViewport)
                gl.viewport(viewport.x, viewport.y, viewport.width, viewport.height)
            },
            get stats(): RendererStats {
                return {
                    renderableCount: scene.count,
                    programCount: ctx.programCache.count,
                    shaderCount: ctx.shaderCache.count,
                    bufferCount: ctx.bufferCount,
                    textureCount: ctx.textureCount,
                    vaoCount: ctx.vaoCount,
                }
            },
            dispose: () => {
                scene.clear()
            }
        }
    }
}

export default Renderer