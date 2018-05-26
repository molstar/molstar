/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { Viewport } from 'mol-view/camera/util';
import { Camera } from 'mol-view/camera/base';

import Scene from './scene';
import { Context } from './webgl/context';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { Renderable } from './renderable';
import { Color } from 'mol-util/color';
import { ValueCell } from 'mol-util';
import { RenderableValues, GlobalUniformValues } from './renderable/schema';
import { RenderObject } from './render-object';

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
    setClearColor: (color: Color) => void

    stats: RendererStats
    dispose: () => void
}

function getPixelRatio() {
    return (typeof window !== 'undefined') ? window.devicePixelRatio : 1
}

export const DefaultRendererProps = {
    clearColor: 0x000000 as Color,
    viewport: Viewport.create(0, 0, 0, 0)
}
export type RendererProps = Partial<typeof DefaultRendererProps>

namespace Renderer {
    export function create(ctx: Context, camera: Camera, props: RendererProps = {}): Renderer {
        const { gl } = ctx
        let { clearColor, viewport: _viewport } = { ...DefaultRendererProps, ...props }
        const scene = Scene.create(ctx)

        const model = Mat4.identity()
        const viewport = Viewport.clone(_viewport)
        const pixelRatio = getPixelRatio()

        // const lightPosition = Vec3.create(0, 0, -100)
        const lightColor = Vec3.create(1.0, 1.0, 1.0)
        const lightAmbient = Vec3.create(0.5, 0.5, 0.5)

        function setClearColor(color: Color) {
            const [ r, g, b ] = Color.toRgbNormalized(color)
            gl.clearColor(r, g, b, 1.0)
        }
        setClearColor(clearColor)

        const globalUniforms: GlobalUniformValues = {
            uModel: ValueCell.create(Mat4.clone(model)),
            uView: ValueCell.create(Mat4.clone(camera.view)),
            uProjection: ValueCell.create(Mat4.clone(camera.projection)),

            uPixelRatio: ValueCell.create(pixelRatio),
            uViewportHeight: ValueCell.create(viewport.height),

            uLightColor: ValueCell.create(Vec3.clone(lightColor)),
            uLightAmbient: ValueCell.create(Vec3.clone(lightAmbient))
        }

        let currentProgramId = -1
        const drawObject = (r: Renderable<RenderableValues>) => {
            if (r.state.visible) {
                if (currentProgramId !== r.program.id) {
                    r.program.use()
                    r.program.setUniforms(globalUniforms)
                    currentProgramId = r.program.id
                }
                if (r.values.dDoubleSided.ref.value) {
                    gl.disable(gl.CULL_FACE)
                } else {
                    gl.enable(gl.CULL_FACE)
                }

                if (r.values.dFlipSided.ref.value) {
                    gl.frontFace(gl.CW)
                    gl.cullFace(gl.FRONT)
                } else {
                    gl.frontFace(gl.CCW)
                    gl.cullFace(gl.BACK)
                }

                gl.depthMask(r.state.depthMask)

                r.draw()
            }
        }

        const draw = () => {
            ValueCell.update(globalUniforms.uView, camera.view)
            ValueCell.update(globalUniforms.uProjection, camera.projection)

            currentProgramId = -1

            gl.depthMask(true)
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

            gl.disable(gl.BLEND)
            gl.enable(gl.DEPTH_TEST)
            scene.eachOpaque(drawObject)

            gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
            gl.enable(gl.BLEND)
            scene.eachTransparent(drawObject)
        }

        return {
            add: (o: RenderObject) => {
                scene.add(o)
            },
            remove: (o: RenderObject) => {
                scene.remove(o)
            },
            update: () => {
                scene.forEach((r, o) => r.update())
            },
            clear: () => {
                scene.clear()
            },
            draw,

            setClearColor,
            setViewport: (newViewport: Viewport) => {
                Viewport.copy(viewport, newViewport)
                gl.viewport(viewport.x, viewport.y, viewport.width, viewport.height)
                ValueCell.update(globalUniforms.uViewportHeight, viewport.height)
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