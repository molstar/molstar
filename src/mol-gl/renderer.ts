/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { Viewport } from 'mol-view/camera/util';
import { Camera } from 'mol-view/camera/base';

import Scene from './scene';
import { Context, createImageData } from './webgl/context';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { Renderable } from './renderable';
import { Color } from 'mol-util/color';
import { ValueCell } from 'mol-util';
import { RenderableValues, GlobalUniformValues } from './renderable/schema';
import { RenderVariant } from './webgl/render-item';

export interface RendererStats {
    programCount: number
    shaderCount: number

    bufferCount: number
    framebufferCount: number
    renderbufferCount: number
    textureCount: number
    vaoCount: number
}

interface Renderer {
    render: (scene: Scene, variant: RenderVariant) => void

    setViewport: (viewport: Viewport) => void
    setClearColor: (color: Color) => void
    getImageData: () => ImageData

    stats: RendererStats
    dispose: () => void
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

        const model = Mat4.identity()
        const viewport = Viewport.clone(_viewport)

        // const lightPosition = Vec3.create(0, 0, -100)
        const lightColor = Vec3.create(1.0, 1.0, 1.0)
        const lightAmbient = Vec3.create(0.5, 0.5, 0.5)
        const highlightColor = Vec3.create(1.0, 0.4, 0.6)
        const selectColor = Vec3.create(0.2, 1.0, 0.1)

        function setClearColor(color: Color) {
            const [ r, g, b ] = Color.toRgbNormalized(color)
            gl.clearColor(r, g, b, 1.0)
        }
        setClearColor(clearColor)

        const globalUniforms: GlobalUniformValues = {
            uModel: ValueCell.create(Mat4.clone(model)),
            uView: ValueCell.create(Mat4.clone(camera.view)),
            uProjection: ValueCell.create(Mat4.clone(camera.projection)),

            uPixelRatio: ValueCell.create(ctx.pixelRatio),
            uViewportHeight: ValueCell.create(viewport.height),

            uLightColor: ValueCell.create(Vec3.clone(lightColor)),
            uLightAmbient: ValueCell.create(Vec3.clone(lightAmbient)),

            uHighlightColor: ValueCell.create(Vec3.clone(highlightColor)),
            uSelectColor: ValueCell.create(Vec3.clone(selectColor))
        }

        let currentProgramId = -1
        const renderObject = (r: Renderable<RenderableValues>, variant: RenderVariant) => {
            const program = r.getProgram(variant)
            if (r.state.visible) {
                if (currentProgramId !== program.id) {
                    program.use()
                    program.setUniforms(globalUniforms)
                    currentProgramId = program.id
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

                r.render(variant)
            }
        }

        const render = (scene: Scene, variant: RenderVariant) => {
            ValueCell.update(globalUniforms.uView, camera.view)
            ValueCell.update(globalUniforms.uProjection, camera.projection)

            currentProgramId = -1

            gl.depthMask(true)
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

            gl.disable(gl.BLEND)
            gl.enable(gl.DEPTH_TEST)
            scene.eachOpaque((r) => renderObject(r, variant))

            gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
            gl.enable(gl.BLEND)
            scene.eachTransparent((r) => renderObject(r, variant))

            gl.finish()
        }

        return {
            render,

            setClearColor,
            setViewport: (newViewport: Viewport) => {
                Viewport.copy(viewport, newViewport)
                gl.viewport(viewport.x, viewport.y, viewport.width, viewport.height)
                ValueCell.update(globalUniforms.uViewportHeight, viewport.height)
            },
            getImageData: () => {
                const { width, height } = viewport
                const buffer = new Uint8Array(width * height * 4)
                ctx.unbindFramebuffer()
                ctx.readPixels(0, 0, width, height, buffer)
                return createImageData(buffer, width, height)
            },

            get stats(): RendererStats {
                console.log(ctx)
                return {
                    programCount: ctx.programCache.count,
                    shaderCount: ctx.shaderCache.count,

                    bufferCount: ctx.bufferCount,
                    framebufferCount: ctx.framebufferCount,
                    renderbufferCount: ctx.renderbufferCount,
                    textureCount: ctx.textureCount,
                    vaoCount: ctx.vaoCount,
                }
            },
            dispose: () => {
                // TODO
            }
        }
    }
}

export default Renderer