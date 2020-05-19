/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Viewport } from '../mol-canvas3d/camera/util';
import { Camera } from '../mol-canvas3d/camera';

import Scene from './scene';
import { WebGLContext } from './webgl/context';
import { Mat4, Vec3, Vec4, Vec2, Quat } from '../mol-math/linear-algebra';
import { Renderable } from './renderable';
import { Color } from '../mol-util/color';
import { ValueCell, deepEqual } from '../mol-util';
import { RenderableValues, GlobalUniformValues, BaseValues } from './renderable/schema';
import { GraphicsRenderVariant } from './webgl/render-item';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { Clipping } from '../mol-theme/clipping';
import { stringToWords } from '../mol-util/string';
import { Transparency } from '../mol-theme/transparency';

export interface RendererStats {
    programCount: number
    shaderCount: number

    attributeCount: number
    elementsCount: number
    framebufferCount: number
    renderbufferCount: number
    textureCount: number
    vertexArrayCount: number

    drawCount: number
    instanceCount: number
    instancedDrawCount: number
}

interface Renderer {
    readonly stats: RendererStats
    readonly props: Readonly<RendererProps>

    clear: (transparentBackground: boolean) => void
    render: (scene: Scene, camera: Camera, variant: GraphicsRenderVariant, clear: boolean, transparentBackground: boolean) => void
    setProps: (props: Partial<RendererProps>) => void
    setViewport: (x: number, y: number, width: number, height: number) => void
    dispose: () => void
}

export const RendererParams = {
    backgroundColor: PD.Color(Color(0x000000), { description: 'Background color of the 3D canvas' }),

    // the following are general 'material' parameters
    pickingAlphaThreshold: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }, { description: 'The minimum opacity value needed for an object to be pickable.' }),
    transparencyVariant: PD.Select('single', PD.arrayToOptions<Transparency.Variant>(['single', 'multi'])),

    interiorDarkening: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }),
    interiorColorFlag: PD.Boolean(true, { label: 'Use Interior Color' }),
    interiorColor: PD.Color(Color.fromNormalizedRgb(0.3, 0.3, 0.3)),

    highlightColor: PD.Color(Color.fromNormalizedRgb(1.0, 0.4, 0.6)),
    selectColor: PD.Color(Color.fromNormalizedRgb(0.2, 1.0, 0.1)),

    style: PD.MappedStatic('matte', {
        custom: PD.Group({
            lightIntensity: PD.Numeric(0.6, { min: 0.0, max: 1.0, step: 0.01 }),
            ambientIntensity: PD.Numeric(0.4, { min: 0.0, max: 1.0, step: 0.01 }),
            metalness: PD.Numeric(0.0, { min: 0.0, max: 1.0, step: 0.01 }),
            roughness: PD.Numeric(1.0, { min: 0.0, max: 1.0, step: 0.01 }),
            reflectivity: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }),
        }, { isExpanded: true }),
        flat: PD.Group({}),
        matte: PD.Group({}),
        glossy: PD.Group({}),
        metallic: PD.Group({}),
        plastic: PD.Group({}),
    }, { label: 'Lighting', description: 'Style in which the 3D scene is rendered/lighted' }),

    clip: PD.Group({
        variant: PD.Select('instance', PD.arrayToOptions<Clipping.Variant>(['instance', 'pixel'])),
        objects: PD.ObjectList({
            type: PD.Select('plane', PD.objectToOptions(Clipping.Type, t => stringToWords(t))),
            position: PD.Vec3(Vec3()),
            rotation: PD.Group({
                axis: PD.Vec3(Vec3.create(1, 0, 0)),
                angle: PD.Numeric(0, { min: -180, max: 180, step: 0.1 }),
            }, { isExpanded: true }),
            scale: PD.Vec3(Vec3.create(1, 1, 1)),
        }, o => stringToWords(o.type))
    })
};
export type RendererProps = PD.Values<typeof RendererParams>

function getStyle(props: RendererProps['style']) {
    switch (props.name) {
        case 'custom':
            return props.params;
        case 'flat':
            return {
                lightIntensity: 0, ambientIntensity: 1,
                metalness: 0, roughness: 0.4, reflectivity: 0.5
            };
        case 'matte':
            return {
                lightIntensity: 0.6, ambientIntensity: 0.4,
                metalness: 0, roughness: 1, reflectivity: 0.5
            };
        case 'glossy':
            return {
                lightIntensity: 0.6, ambientIntensity: 0.4,
                metalness: 0, roughness: 0.4, reflectivity: 0.5
            };
        case 'metallic':
            return {
                lightIntensity: 0.6, ambientIntensity: 0.4,
                metalness: 0.4, roughness: 0.6, reflectivity: 0.5
            };
        case 'plastic':
            return {
                lightIntensity: 0.6, ambientIntensity: 0.4,
                metalness: 0, roughness: 0.2, reflectivity: 0.5
            };
    }
}

type Clip = {
    variant: Clipping.Variant
    objects: {
        count: number
        type: number[]
        position: number[]
        rotation: number[]
        scale: number[]
    }
}

const tmpQuat = Quat();
function getClip(props: RendererProps['clip'], clip?: Clip): Clip {
    const { type, position, rotation, scale } = clip?.objects || {
        type: (new Array(5)).fill(1),
        position: (new Array(5 * 3)).fill(0),
        rotation: (new Array(5 * 4)).fill(0),
        scale: (new Array(5 * 3)).fill(1),
    };
    for (let i = 0, il = props.objects.length; i < il; ++i) {
        const p = props.objects[i];
        type[i] = Clipping.Type[p.type];
        Vec3.toArray(p.position, position, i * 3);
        Quat.toArray(Quat.setAxisAngle(tmpQuat, p.rotation.axis, p.rotation.angle), rotation, i * 4);
        Vec3.toArray(p.scale, scale, i * 3);
    }
    return {
        variant: props.variant,
        objects: { count: props.objects.length, type, position, rotation, scale }
    };
}

namespace Renderer {
    export function create(ctx: WebGLContext, props: Partial<RendererProps> = {}): Renderer {
        const { gl, state, stats } = ctx;
        const p = PD.merge(RendererParams, PD.getDefaultValues(RendererParams), props);
        const style = getStyle(p.style);
        const clip = getClip(p.clip);

        const viewport = Viewport();
        const bgColor = Color.toVec3Normalized(Vec3(), p.backgroundColor);

        const view = Mat4();
        const invView = Mat4();
        const modelView = Mat4();
        const invModelView = Mat4();
        const invProjection = Mat4();
        const modelViewProjection = Mat4();
        const invModelViewProjection = Mat4();

        const viewOffset = Vec2();

        const globalUniforms: GlobalUniformValues = {
            uModel: ValueCell.create(Mat4.identity()),
            uView: ValueCell.create(view),
            uInvView: ValueCell.create(invView),
            uModelView: ValueCell.create(modelView),
            uInvModelView: ValueCell.create(invModelView),
            uInvProjection: ValueCell.create(invProjection),
            uProjection: ValueCell.create(Mat4()),
            uModelViewProjection: ValueCell.create(modelViewProjection),
            uInvModelViewProjection: ValueCell.create(invModelViewProjection),

            uIsOrtho: ValueCell.create(1),
            uViewOffset: ValueCell.create(viewOffset),

            uPixelRatio: ValueCell.create(ctx.pixelRatio),
            uViewportHeight: ValueCell.create(viewport.height),
            uViewport: ValueCell.create(Viewport.toVec4(Vec4(), viewport)),

            uCameraPosition: ValueCell.create(Vec3()),
            uNear: ValueCell.create(1),
            uFar: ValueCell.create(10000),
            uFogNear: ValueCell.create(1),
            uFogFar: ValueCell.create(10000),
            uFogColor: ValueCell.create(bgColor),
            uTransparentBackground: ValueCell.create(0),

            uClipObjectType: ValueCell.create(clip.objects.type),
            uClipObjectPosition: ValueCell.create(clip.objects.position),
            uClipObjectRotation: ValueCell.create(clip.objects.rotation),
            uClipObjectScale: ValueCell.create(clip.objects.scale),

            // the following are general 'material' uniforms
            uLightIntensity: ValueCell.create(style.lightIntensity),
            uAmbientIntensity: ValueCell.create(style.ambientIntensity),

            uMetalness: ValueCell.create(style.metalness),
            uRoughness: ValueCell.create(style.roughness),
            uReflectivity: ValueCell.create(style.reflectivity),

            uPickingAlphaThreshold: ValueCell.create(p.pickingAlphaThreshold),

            uInteriorDarkening: ValueCell.create(p.interiorDarkening),
            uInteriorColorFlag: ValueCell.create(p.interiorColorFlag ? 1 : 0),
            uInteriorColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.interiorColor)),

            uHighlightColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.highlightColor)),
            uSelectColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.selectColor)),
        };
        const globalUniformList = Object.entries(globalUniforms);

        let globalUniformsNeedUpdate = true;

        const renderObject = (r: Renderable<RenderableValues & BaseValues>, variant: GraphicsRenderVariant) => {
            if (!r.state.visible || (!r.state.pickable && variant[0] === 'p')) {
                return;
            }

            let definesNeedUpdate = false;
            if (r.values.dClipObjectCount.ref.value !== clip.objects.count) {
                ValueCell.update(r.values.dClipObjectCount, clip.objects.count);
                definesNeedUpdate = true;
            }
            if (r.values.dClipVariant.ref.value !== clip.variant) {
                ValueCell.update(r.values.dClipVariant, clip.variant);
                definesNeedUpdate = true;
            }
            if (r.values.dTransparencyVariant.ref.value !== p.transparencyVariant) {
                ValueCell.update(r.values.dTransparencyVariant, p.transparencyVariant);
                definesNeedUpdate = true;
            }
            if (definesNeedUpdate) r.update();

            const program = r.getProgram(variant);
            if (state.currentProgramId !== program.id) {
                // console.log('new program')
                globalUniformsNeedUpdate = true;
                program.use();
            }

            if (globalUniformsNeedUpdate) {
                // console.log('globalUniformsNeedUpdate')
                program.setUniforms(globalUniformList);
                globalUniformsNeedUpdate = false;
            }

            if (r.values.dDoubleSided) {
                if (r.values.dDoubleSided.ref.value) {
                    state.disable(gl.CULL_FACE);
                } else {
                    state.enable(gl.CULL_FACE);
                }
            } else {
                // webgl default
                state.disable(gl.CULL_FACE);
            }

            if (r.values.dFlipSided) {
                if (r.values.dFlipSided.ref.value) {
                    state.frontFace(gl.CW);
                    state.cullFace(gl.FRONT);
                } else {
                    state.frontFace(gl.CCW);
                    state.cullFace(gl.BACK);
                }
            } else {
                // webgl default
                state.frontFace(gl.CCW);
                state.cullFace(gl.BACK);
            }

            if (variant === 'color') {
                state.depthMask(r.state.writeDepth);
            }

            r.render(variant);
        };

        const render = (scene: Scene, camera: Camera, variant: GraphicsRenderVariant, clear: boolean, transparentBackground: boolean) => {
            ValueCell.update(globalUniforms.uModel, scene.view);
            ValueCell.update(globalUniforms.uView, camera.view);
            ValueCell.update(globalUniforms.uInvView, Mat4.invert(invView, camera.view));
            ValueCell.update(globalUniforms.uModelView, Mat4.mul(modelView, scene.view, camera.view));
            ValueCell.update(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView));
            ValueCell.update(globalUniforms.uProjection, camera.projection);
            ValueCell.update(globalUniforms.uInvProjection, Mat4.invert(invProjection, camera.projection));
            ValueCell.update(globalUniforms.uModelViewProjection, Mat4.mul(modelViewProjection, modelView, camera.projection));
            ValueCell.update(globalUniforms.uInvModelViewProjection, Mat4.invert(invModelViewProjection, modelViewProjection));

            ValueCell.update(globalUniforms.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
            ValueCell.update(globalUniforms.uViewOffset, camera.viewOffset.enabled ? Vec2.set(viewOffset, camera.viewOffset.offsetX * 16, camera.viewOffset.offsetY * 16) : Vec2.set(viewOffset, 0, 0));

            ValueCell.update(globalUniforms.uCameraPosition, camera.state.position);
            ValueCell.update(globalUniforms.uFar, camera.far);
            ValueCell.update(globalUniforms.uNear, camera.near);
            ValueCell.update(globalUniforms.uFogFar, camera.fogFar);
            ValueCell.update(globalUniforms.uFogNear, camera.fogNear);

            ValueCell.update(globalUniforms.uTransparentBackground, transparentBackground ? 1 : 0);

            globalUniformsNeedUpdate = true;
            state.currentRenderItemId = -1;

            const { renderables } = scene;

            state.disable(gl.SCISSOR_TEST);
            state.disable(gl.BLEND);
            state.colorMask(true, true, true, true);
            state.enable(gl.DEPTH_TEST);

            if (clear) {
                state.depthMask(true);
                if (variant === 'color') {
                    state.clearColor(bgColor[0], bgColor[1], bgColor[2], transparentBackground ? 0 : 1);
                } else {
                    state.clearColor(1, 1, 1, 1);
                }
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
            }

            if (variant === 'color') {
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    const r = renderables[i];
                    if (r.state.opaque) renderObject(r, variant);
                }

                state.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE);
                state.enable(gl.BLEND);
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    const r = renderables[i];
                    if (!r.state.opaque && r.state.writeDepth) renderObject(r, variant);
                }
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    const r = renderables[i];
                    if (!r.state.opaque && !r.state.writeDepth) renderObject(r, variant);
                }
            } else { // picking & depth
                for (let i = 0, il = renderables.length; i < il; ++i) {
                    renderObject(renderables[i], variant);
                }
            }

            gl.finish();
        };

        return {
            clear: (transparentBackground: boolean) => {
                state.depthMask(true);
                state.colorMask(true, true, true, true);
                state.clearColor(bgColor[0], bgColor[1], bgColor[2], transparentBackground ? 0 : 1);
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
            },
            render,

            setProps: (props: Partial<RendererProps>) => {
                if (props.backgroundColor !== undefined && props.backgroundColor !== p.backgroundColor) {
                    p.backgroundColor = props.backgroundColor;
                    Color.toVec3Normalized(bgColor, p.backgroundColor);
                    ValueCell.update(globalUniforms.uFogColor, Vec3.copy(globalUniforms.uFogColor.ref.value, bgColor));
                }

                if (props.pickingAlphaThreshold !== undefined && props.pickingAlphaThreshold !== p.pickingAlphaThreshold) {
                    p.pickingAlphaThreshold = props.pickingAlphaThreshold;
                    ValueCell.update(globalUniforms.uPickingAlphaThreshold, p.pickingAlphaThreshold);
                }
                if (props.transparencyVariant !== undefined && props.transparencyVariant !== p.transparencyVariant) {
                    p.transparencyVariant = props.transparencyVariant;
                }

                if (props.interiorDarkening !== undefined && props.interiorDarkening !== p.interiorDarkening) {
                    p.interiorDarkening = props.interiorDarkening;
                    ValueCell.update(globalUniforms.uInteriorDarkening, p.interiorDarkening);
                }
                if (props.interiorColorFlag !== undefined && props.interiorColorFlag !== p.interiorColorFlag) {
                    p.interiorColorFlag = props.interiorColorFlag;
                    ValueCell.update(globalUniforms.uInteriorColorFlag, p.interiorColorFlag ? 1 : 0);
                }
                if (props.interiorColor !== undefined && props.interiorColor !== p.interiorColor) {
                    p.interiorColor = props.interiorColor;
                    ValueCell.update(globalUniforms.uInteriorColor, Color.toVec3Normalized(globalUniforms.uInteriorColor.ref.value, p.interiorColor));
                }

                if (props.highlightColor !== undefined && props.highlightColor !== p.highlightColor) {
                    p.highlightColor = props.highlightColor;
                    ValueCell.update(globalUniforms.uHighlightColor, Color.toVec3Normalized(globalUniforms.uHighlightColor.ref.value, p.highlightColor));
                }
                if (props.selectColor !== undefined && props.selectColor !== p.selectColor) {
                    p.selectColor = props.selectColor;
                    ValueCell.update(globalUniforms.uSelectColor, Color.toVec3Normalized(globalUniforms.uSelectColor.ref.value, p.selectColor));
                }

                if (props.style !== undefined) {
                    p.style = props.style;
                    Object.assign(style, getStyle(props.style));
                    ValueCell.updateIfChanged(globalUniforms.uLightIntensity, style.lightIntensity);
                    ValueCell.updateIfChanged(globalUniforms.uAmbientIntensity, style.ambientIntensity);
                    ValueCell.updateIfChanged(globalUniforms.uMetalness, style.metalness);
                    ValueCell.updateIfChanged(globalUniforms.uRoughness, style.roughness);
                    ValueCell.updateIfChanged(globalUniforms.uReflectivity, style.reflectivity);
                }

                if (props.clip !== undefined && !deepEqual(props.clip, p.clip)) {
                    p.clip = props.clip;
                    Object.assign(clip, getClip(props.clip, clip));
                    ValueCell.update(globalUniforms.uClipObjectPosition, clip.objects.position);
                    ValueCell.update(globalUniforms.uClipObjectRotation, clip.objects.rotation);
                    ValueCell.update(globalUniforms.uClipObjectScale, clip.objects.scale);
                    ValueCell.update(globalUniforms.uClipObjectType, clip.objects.type);
                }
            },
            setViewport: (x: number, y: number, width: number, height: number) => {
                gl.viewport(x, y, width, height);
                if (x !== viewport.x || y !== viewport.y || width !== viewport.width || height !== viewport.height) {
                    Viewport.set(viewport, x, y, width, height);
                    ValueCell.update(globalUniforms.uViewportHeight, height);
                    ValueCell.update(globalUniforms.uViewport, Vec4.set(globalUniforms.uViewport.ref.value, x, y, width, height));
                }
            },

            get props() {
                return p;
            },
            get stats(): RendererStats {
                return {
                    programCount: ctx.stats.resourceCounts.program,
                    shaderCount: ctx.stats.resourceCounts.shader,

                    attributeCount: ctx.stats.resourceCounts.attribute,
                    elementsCount: ctx.stats.resourceCounts.elements,
                    framebufferCount: ctx.stats.resourceCounts.framebuffer,
                    renderbufferCount: ctx.stats.resourceCounts.renderbuffer,
                    textureCount: ctx.stats.resourceCounts.texture,
                    vertexArrayCount: ctx.stats.resourceCounts.vertexArray,

                    drawCount: stats.drawCount,
                    instanceCount: stats.instanceCount,
                    instancedDrawCount: stats.instancedDrawCount,
                };
            },
            dispose: () => {
                // TODO
            }
        };
    }
}

export default Renderer;