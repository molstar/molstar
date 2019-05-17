/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Viewport } from 'mol-canvas3d/camera/util';
import { Camera } from 'mol-canvas3d/camera';

import Scene from './scene';
import { WebGLContext } from './webgl/context';
import { Mat4, Vec3, Vec4 } from 'mol-math/linear-algebra';
import { Renderable } from './renderable';
import { Color } from 'mol-util/color';
import { ValueCell } from 'mol-util';
import { RenderableValues, GlobalUniformValues, BaseValues } from './renderable/schema';
import { GraphicsRenderVariant } from './webgl/render-item';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { deepClone } from 'mol-util/object';

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
    readonly props: Readonly<RendererProps>

    clear: () => void
    render: (scene: Scene, variant: GraphicsRenderVariant, clear: boolean) => void
    setProps: (props: Partial<RendererProps>) => void
    setViewport: (x: number, y: number, width: number, height: number) => void
    dispose: () => void
}

export const RendererParams = {
    backgroundColor: PD.Color(Color(0x000000)),
    pickingAlphaThreshold: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }, { description: 'The minimum opacity value needed for an object to be pickable.' }),

    lightIntensity: PD.Numeric(0.6, { min: 0.0, max: 1.0, step: 0.01 }),
    ambientIntensity: PD.Numeric(0.4, { min: 0.0, max: 1.0, step: 0.01 }),

    metalness: PD.Numeric(0.0, { min: 0.0, max: 1.0, step: 0.01 }),
    roughness: PD.Numeric(0.4, { min: 0.0, max: 1.0, step: 0.01 }),
    reflectivity: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }),
}
export type RendererProps = PD.Values<typeof RendererParams>

namespace Renderer {
    export function create(ctx: WebGLContext, camera: Camera, props: Partial<RendererProps> = {}): Renderer {
        const { gl, state, stats } = ctx
        const p = deepClone({ ...PD.getDefaultValues(RendererParams), ...props })

        const viewport = Viewport()
        const bgColor = Color.toVec3Normalized(Vec3(), p.backgroundColor)

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

            uIsOrtho: ValueCell.create(camera.state.mode === 'orthographic' ? 1 : 0),
            uPixelRatio: ValueCell.create(ctx.pixelRatio),
            uViewportHeight: ValueCell.create(viewport.height),
            uViewport: ValueCell.create(Viewport.toVec4(Vec4(), viewport)),

            uLightIntensity: ValueCell.create(p.lightIntensity),
            uAmbientIntensity: ValueCell.create(p.ambientIntensity),

            uMetalness: ValueCell.create(p.metalness),
            uRoughness: ValueCell.create(p.roughness),
            uReflectivity: ValueCell.create(p.reflectivity),

            uCameraPosition: ValueCell.create(Vec3.clone(camera.state.position)),
            uNear: ValueCell.create(camera.state.near),
            uFar: ValueCell.create(camera.state.far),
            uFogNear: ValueCell.create(camera.state.fogNear),
            uFogFar: ValueCell.create(camera.state.fogFar),
            uFogColor: ValueCell.create(bgColor),

            uPickingAlphaThreshold: ValueCell.create(p.pickingAlphaThreshold),
        }
        const globalUniformList = Object.entries(globalUniforms)

        let globalUniformsNeedUpdate = true

        const renderObject = (r: Renderable<RenderableValues & BaseValues>, variant: GraphicsRenderVariant) => {
            const program = r.getProgram(variant)
            if (r.state.visible) {
                if (state.currentProgramId !== program.id) {
                    // console.log('new program')
                    globalUniformsNeedUpdate = true
                    program.use()
                }

                if (globalUniformsNeedUpdate) {
                    // console.log('globalUniformsNeedUpdate')
                    program.setUniforms(globalUniformList)
                    globalUniformsNeedUpdate = false
                }

                if (r.values.dDoubleSided) {
                    if (r.values.dDoubleSided.ref.value) {
                        state.disable(gl.CULL_FACE)
                    } else {
                        state.enable(gl.CULL_FACE)
                    }
                } else {
                    // webgl default
                    state.disable(gl.CULL_FACE)
                }

                if (r.values.dFlipSided) {
                    if (r.values.dFlipSided.ref.value) {
                        state.frontFace(gl.CW)
                        state.cullFace(gl.FRONT)
                    } else {
                        state.frontFace(gl.CCW)
                        state.cullFace(gl.BACK)
                    }
                } else {
                    // webgl default
                    state.frontFace(gl.CCW)
                    state.cullFace(gl.BACK)
                }

                r.render(variant)
            }
        }

        const render = (scene: Scene, variant: GraphicsRenderVariant, clear: boolean) => {
            ValueCell.update(globalUniforms.uModel, scene.view)
            ValueCell.update(globalUniforms.uView, camera.view)
            ValueCell.update(globalUniforms.uInvView, Mat4.invert(invView, camera.view))
            ValueCell.update(globalUniforms.uModelView, Mat4.mul(modelView, scene.view, camera.view))
            ValueCell.update(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView))
            ValueCell.update(globalUniforms.uProjection, camera.projection)
            ValueCell.update(globalUniforms.uInvProjection, Mat4.invert(invProjection, camera.projection))
            ValueCell.update(globalUniforms.uModelViewProjection, Mat4.mul(modelViewProjection, modelView, camera.projection))
            ValueCell.update(globalUniforms.uInvModelViewProjection, Mat4.invert(invModelViewProjection, modelViewProjection))

            ValueCell.update(globalUniforms.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0)

            ValueCell.update(globalUniforms.uCameraPosition, camera.state.position)
            ValueCell.update(globalUniforms.uFar, camera.state.far)
            ValueCell.update(globalUniforms.uNear, camera.state.near)
            ValueCell.update(globalUniforms.uFogFar, camera.state.fogFar)
            ValueCell.update(globalUniforms.uFogNear, camera.state.fogNear)

            globalUniformsNeedUpdate = true
            state.currentRenderItemId = -1

            const { renderables } = scene

            state.disable(gl.SCISSOR_TEST)
            state.disable(gl.BLEND)
            state.depthMask(true)
            state.colorMask(true, true, true, true)
            state.enable(gl.DEPTH_TEST)

            if (clear) {
                if (variant === 'color') {
                    state.clearColor(bgColor[0], bgColor[1], bgColor[2], 1.0)
                } else {
                    state.clearColor(1, 1, 1, 1)
                }
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
            }

            if (variant === 'color') {
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    const r = renderables[i]
                    if (r.state.opaque) renderObject(r, variant)
                }

                state.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
                state.enable(gl.BLEND)
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    const r = renderables[i]
                    state.depthMask(r.values.uAlpha.ref.value === 1.0)
                    if (!r.state.opaque) renderObject(r, variant)
                }
            } else { // picking & depth
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    renderObject(renderables[i], variant)
                }
            }

            gl.finish()
        }

        return {
            clear: () => {
                state.depthMask(true)
                state.colorMask(true, true, true, true)
                state.clearColor(bgColor[0], bgColor[1], bgColor[2], 1.0)
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
            },
            render,

            setProps: (props: Partial<RendererProps>) => {
                if (props.pickingAlphaThreshold !== undefined && props.pickingAlphaThreshold !== p.pickingAlphaThreshold) {
                    p.pickingAlphaThreshold = props.pickingAlphaThreshold
                    ValueCell.update(globalUniforms.uPickingAlphaThreshold, p.pickingAlphaThreshold)
                }
                if (props.backgroundColor !== undefined && props.backgroundColor !== p.backgroundColor) {
                    p.backgroundColor = props.backgroundColor
                    Color.toVec3Normalized(bgColor, p.backgroundColor)
                    ValueCell.update(globalUniforms.uFogColor, Vec3.copy(globalUniforms.uFogColor.ref.value, bgColor))
                }
                if (props.lightIntensity !== undefined && props.lightIntensity !== p.lightIntensity) {
                    p.lightIntensity = props.lightIntensity
                    ValueCell.update(globalUniforms.uLightIntensity, p.lightIntensity)
                }
                if (props.ambientIntensity !== undefined && props.ambientIntensity !== p.ambientIntensity) {
                    p.ambientIntensity = props.ambientIntensity
                    ValueCell.update(globalUniforms.uAmbientIntensity, p.ambientIntensity)
                }

                if (props.metalness !== undefined && props.metalness !== p.metalness) {
                    p.metalness = props.metalness
                    ValueCell.update(globalUniforms.uMetalness, p.metalness)
                }
                if (props.roughness !== undefined && props.roughness !== p.roughness) {
                    p.roughness = props.roughness
                    ValueCell.update(globalUniforms.uRoughness, p.roughness)
                }
                if (props.reflectivity !== undefined && props.reflectivity !== p.reflectivity) {
                    p.reflectivity = props.reflectivity
                    ValueCell.update(globalUniforms.uReflectivity, p.reflectivity)
                }
            },
            setViewport: (x: number, y: number, width: number, height: number) => {
                gl.viewport(x, y, width, height)
                if (x !== viewport.x || y !== viewport.y || width !== viewport.width || height !== viewport.height) {
                    Viewport.set(viewport, x, y, width, height)
                    ValueCell.update(globalUniforms.uViewportHeight, height)
                    ValueCell.update(globalUniforms.uViewport, Vec4.set(globalUniforms.uViewport.ref.value, x, y, width, height))
                }
            },

            get props() {
                return p
            },
            get stats(): RendererStats {
                return {
                    programCount: ctx.programCache.count,
                    shaderCount: ctx.shaderCache.count,

                    bufferCount: stats.bufferCount,
                    framebufferCount: stats.framebufferCount,
                    renderbufferCount: stats.renderbufferCount,
                    textureCount: stats.textureCount,
                    vaoCount: stats.vaoCount,

                    drawCount: stats.drawCount,
                    instanceCount: stats.instanceCount,
                    instancedDrawCount: stats.instancedDrawCount,
                }
            },
            dispose: () => {
                // TODO
            }
        }
    }
}

export default Renderer