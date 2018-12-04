/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { Viewport } from 'mol-canvas3d/camera/util';
import { Camera } from 'mol-canvas3d/camera';

import Scene from './scene';
import { WebGLContext, createImageData } from './webgl/context';
import { Mat4, Vec3, Vec4 } from 'mol-math/linear-algebra';
import { Renderable } from './renderable';
import { Color } from 'mol-util/color';
import { ValueCell } from 'mol-util';
import { RenderableValues, GlobalUniformValues, BaseValues } from './renderable/schema';
import { RenderVariant } from './webgl/render-item';

export interface RendererStats {
    programCount: number
    shaderCount: number

    bufferCount: number
    framebufferCount: number
    renderbufferCount: number
    textureCount: number
    vaoCount: number

    drawCount: number
    instanceCount: number
    instancedDrawCount: number
}

interface Renderer {
    readonly stats: RendererStats
    readonly props: RendererProps

    clear: () => void
    render: (scene: Scene, variant: RenderVariant) => void
    setViewport: (x: number, y: number, width: number, height: number) => void
    setClearColor: (color: Color) => void
    setPickingAlphaThreshold: (value: number) => void
    getImageData: () => ImageData
    dispose: () => void
}

export const DefaultRendererProps = {
    clearColor: Color(0x000000),
    viewport: Viewport.create(0, 0, 0, 0),
    pickingAlphaThreshold: 0.5,
}
export type RendererProps = typeof DefaultRendererProps

namespace Renderer {
    export function create(ctx: WebGLContext, camera: Camera, props: Partial<RendererProps> = {}): Renderer {
        const { gl } = ctx
        let { clearColor, viewport: _viewport, pickingAlphaThreshold } = { ...DefaultRendererProps, ...props }

        const viewport = Viewport.clone(_viewport)
        const viewportVec4 = Viewport.toVec4(Vec4.zero(), viewport)

        // const lightPosition = Vec3.create(0, 0, -100)
        const lightColor = Vec3.create(1.0, 1.0, 1.0)
        const lightAmbient = Vec3.create(0.5, 0.5, 0.5)
        const fogColor = Vec3.create(0.0, 0.0, 0.0)

        function setClearColor(color: Color) {
            clearColor = color
            const [ r, g, b ] = Color.toRgbNormalized(color)
            gl.clearColor(r, g, b, 1.0)
            Vec3.set(fogColor, r, g, b)
        }
        setClearColor(clearColor)

        const view = Mat4.clone(camera.view)
        const invView = Mat4.invert(Mat4.identity(), view)
        const modelView = Mat4.clone(camera.view)
        const invModelView = Mat4.invert(Mat4.identity(), modelView)
        const invProjection = Mat4.invert(Mat4.identity(), camera.projection)
        const modelViewProjection = Mat4.mul(Mat4.identity(), modelView, camera.projection)
        const invModelViewProjection = Mat4.invert(Mat4.identity(), modelViewProjection)

        const globalUniforms: GlobalUniformValues = {
            uModel: ValueCell.create(Mat4.identity()),
            uView: ValueCell.create(camera.view),
            uInvView: ValueCell.create(invView),
            uModelView: ValueCell.create(modelView),
            uInvModelView: ValueCell.create(invModelView),
            uInvProjection: ValueCell.create(invProjection),
            uProjection: ValueCell.create(Mat4.clone(camera.projection)),
            uModelViewProjection: ValueCell.create(modelViewProjection),
            uInvModelViewProjection: ValueCell.create(invModelViewProjection),

            uPixelRatio: ValueCell.create(ctx.pixelRatio),
            uViewportHeight: ValueCell.create(viewport.height),
            uViewport: ValueCell.create(viewportVec4),

            uLightColor: ValueCell.create(lightColor),
            uLightAmbient: ValueCell.create(lightAmbient),

            uFogNear: ValueCell.create(camera.state.fogNear),
            uFogFar: ValueCell.create(camera.state.fogFar),
            uFogColor: ValueCell.create(fogColor),

            uPickingAlphaThreshold: ValueCell.create(pickingAlphaThreshold),
        }

        let currentProgramId = -1
        const renderObject = (r: Renderable<RenderableValues & BaseValues>, variant: RenderVariant, opaque: boolean) => {
            if (r.state.opaque !== opaque) return
            const program = r.getProgram(variant)
            if (r.state.visible) {
                if (currentProgramId !== program.id) {
                    program.use()
                    program.setUniforms(globalUniforms)
                    currentProgramId = program.id
                }

                if (r.values.dDoubleSided) {
                    if (r.values.dDoubleSided.ref.value) {
                        gl.disable(gl.CULL_FACE)
                    } else {
                        gl.enable(gl.CULL_FACE)
                    }
                } else {
                    // webgl default
                    gl.disable(gl.CULL_FACE)
                }

                if (r.values.dFlipSided) {
                    if (r.values.dFlipSided.ref.value) {
                        gl.frontFace(gl.CW)
                        gl.cullFace(gl.FRONT)
                    } else {
                        gl.frontFace(gl.CCW)
                        gl.cullFace(gl.BACK)
                    }
                } else {
                    // webgl default
                    gl.frontFace(gl.CCW)
                    gl.cullFace(gl.BACK)
                }

                gl.depthMask(r.state.opaque)

                r.render(variant)
            }
        }

        const render = (scene: Scene, variant: RenderVariant) => {
            ValueCell.update(globalUniforms.uModel, scene.view)
            ValueCell.update(globalUniforms.uView, camera.view)
            ValueCell.update(globalUniforms.uInvView, Mat4.invert(invView, camera.view))
            ValueCell.update(globalUniforms.uModelView, Mat4.mul(modelView, scene.view, camera.view))
            ValueCell.update(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView))
            ValueCell.update(globalUniforms.uProjection, camera.projection)
            ValueCell.update(globalUniforms.uInvProjection, Mat4.invert(invProjection, camera.projection))
            ValueCell.update(globalUniforms.uModelViewProjection, Mat4.mul(modelViewProjection, modelView, camera.projection))
            ValueCell.update(globalUniforms.uInvModelViewProjection, Mat4.invert(invModelViewProjection, modelViewProjection))

            ValueCell.update(globalUniforms.uFogFar, camera.state.fogFar)
            ValueCell.update(globalUniforms.uFogNear, camera.state.fogNear)

            const { renderables } = scene
            currentProgramId = -1

            gl.disable(gl.BLEND)
            gl.enable(gl.DEPTH_TEST)
            for (let i = 0, il = renderables.length; i < il; ++i) {
                renderObject(renderables[i], variant, true)
            }

            gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
            gl.enable(gl.BLEND)
            for (let i = 0, il = renderables.length; i < il; ++i) {
                renderObject(renderables[i], variant, false)
            }

            gl.finish()
        }

        return {
            clear: () => {
                gl.depthMask(true)
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
            },
            render,

            setClearColor,
            setPickingAlphaThreshold: (value: number) => {
                pickingAlphaThreshold = value
                ValueCell.update(globalUniforms.uPickingAlphaThreshold, pickingAlphaThreshold)
            },
            setViewport: (x: number, y: number, width: number, height: number) => {
                Viewport.set(viewport, x, y, width, height)
                gl.viewport(x, y, width, height)
                ValueCell.update(globalUniforms.uViewportHeight, height)
                ValueCell.update(globalUniforms.uViewport, Vec4.set(viewportVec4, x, y, width, height))
            },
            getImageData: () => {
                const { width, height } = viewport
                const buffer = new Uint8Array(width * height * 4)
                ctx.unbindFramebuffer()
                ctx.readPixels(0, 0, width, height, buffer)
                return createImageData(buffer, width, height)
            },

            get props() {
                return {
                    clearColor,
                    pickingAlphaThreshold,
                    viewport
                }
            },
            get stats(): RendererStats {
                return {
                    programCount: ctx.programCache.count,
                    shaderCount: ctx.shaderCache.count,

                    bufferCount: ctx.bufferCount,
                    framebufferCount: ctx.framebufferCount,
                    renderbufferCount: ctx.renderbufferCount,
                    textureCount: ctx.textureCount,
                    vaoCount: ctx.vaoCount,

                    drawCount: ctx.drawCount,
                    instanceCount: ctx.instanceCount,
                    instancedDrawCount: ctx.instancedDrawCount,
                }
            },
            dispose: () => {
                // TODO
            }
        }
    }
}

export default Renderer