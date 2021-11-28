/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Viewport } from '../mol-canvas3d/camera/util';
import { ICamera } from '../mol-canvas3d/camera';
import { Scene } from './scene';
import { WebGLContext } from './webgl/context';
import { Mat4, Vec3, Vec4, Vec2, Quat } from '../mol-math/linear-algebra';
import { GraphicsRenderable } from './renderable';
import { Color } from '../mol-util/color';
import { ValueCell, deepEqual } from '../mol-util';
import { GlobalUniformValues } from './renderable/schema';
import { GraphicsRenderVariant } from './webgl/render-item';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { Clipping } from '../mol-theme/clipping';
import { stringToWords } from '../mol-util/string';
import { degToRad } from '../mol-math/misc';
import { createNullTexture, Texture, Textures } from './webgl/texture';
import { arrayMapUpsert } from '../mol-util/array';
import { clamp } from '../mol-math/interpolate';

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

    clear: (toBackgroundColor: boolean) => void
    clearDepth: () => void
    update: (camera: ICamera) => void

    renderPick: (group: Scene.Group, camera: ICamera, variant: GraphicsRenderVariant, depthTexture: Texture | null) => void
    renderDepth: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderMarkingDepth: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderMarkingMask: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlended: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlendedOpaque: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlendedTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlendedVolumeOpaque: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlendedVolumeTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderWboitOpaque: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderWboitTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void

    setProps: (props: Partial<RendererProps>) => void
    setViewport: (x: number, y: number, width: number, height: number) => void
    setTransparentBackground: (value: boolean) => void
    setDrawingBufferSize: (width: number, height: number) => void
    setPixelRatio: (value: number) => void

    dispose: () => void
}

export const RendererParams = {
    backgroundColor: PD.Color(Color(0x000000), { description: 'Background color of the 3D canvas' }),

    pickingAlphaThreshold: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }, { description: 'The minimum opacity value needed for an object to be pickable.' }),

    interiorDarkening: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }),
    interiorColorFlag: PD.Boolean(true, { label: 'Use Interior Color' }),
    interiorColor: PD.Color(Color.fromNormalizedRgb(0.3, 0.3, 0.3)),

    highlightColor: PD.Color(Color.fromNormalizedRgb(1.0, 0.4, 0.6)),
    selectColor: PD.Color(Color.fromNormalizedRgb(0.2, 1.0, 0.1)),
    highlightStrength: PD.Numeric(0.7, { min: 0.0, max: 1.0, step: 0.1 }),
    selectStrength: PD.Numeric(0.7, { min: 0.0, max: 1.0, step: 0.1 }),
    markerPriority: PD.Select(1, [[1, 'Highlight'], [2, 'Select']]),

    xrayEdgeFalloff: PD.Numeric(1, { min: 0.0, max: 3.0, step: 0.1 }),

    light: PD.ObjectList({
        inclination: PD.Numeric(180, { min: 0, max: 180, step: 1 }),
        azimuth: PD.Numeric(0, { min: 0, max: 360, step: 1 }),
        color: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
        intensity: PD.Numeric(0.6, { min: 0.0, max: 1.0, step: 0.01 }),
    }, o => Color.toHexString(o.color), { defaultValue: [{
        inclination: 180,
        azimuth: 0,
        color: Color.fromNormalizedRgb(1.0, 1.0, 1.0),
        intensity: 0.6
    }] }),
    ambientColor: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
    ambientIntensity: PD.Numeric(0.4, { min: 0.0, max: 1.0, step: 0.01 }),

    clip: PD.Group({
        variant: PD.Select('instance', PD.arrayToOptions<Clipping.Variant>(['instance', 'pixel'])),
        objects: PD.ObjectList({
            type: PD.Select('plane', PD.objectToOptions(Clipping.Type, t => stringToWords(t))),
            invert: PD.Boolean(false),
            position: PD.Vec3(Vec3()),
            rotation: PD.Group({
                axis: PD.Vec3(Vec3.create(1, 0, 0)),
                angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { description: 'Angle in Degrees' }),
            }, { isExpanded: true }),
            scale: PD.Vec3(Vec3.create(1, 1, 1)),
        }, o => stringToWords(o.type))
    })
};
export type RendererProps = PD.Values<typeof RendererParams>

type Light = {
    count: number
    direction: number[]
    color: number[]
}

const tmpDir = Vec3();
const tmpColor = Vec3();
function getLight(props: RendererProps['light'], light?: Light): Light {
    const { direction, color } = light || {
        direction: (new Array(5 * 3)).fill(0),
        color: (new Array(5 * 3)).fill(0),
    };
    for (let i = 0, il = props.length; i < il; ++i) {
        const p = props[i];
        Vec3.directionFromSpherical(tmpDir, degToRad(p.inclination), degToRad(p.azimuth), 1);
        Vec3.toArray(tmpDir, direction, i * 3);
        Vec3.scale(tmpColor, Color.toVec3Normalized(tmpColor, p.color), p.intensity);
        Vec3.toArray(tmpColor, color, i * 3);
    }
    return { count: props.length, direction, color };
}

type Clip = {
    variant: Clipping.Variant
    objects: {
        count: number
        type: number[]
        invert: boolean[]
        position: number[]
        rotation: number[]
        scale: number[]
    }
}

const tmpQuat = Quat();
function getClip(props: RendererProps['clip'], clip?: Clip): Clip {
    const { type, invert, position, rotation, scale } = clip?.objects || {
        type: (new Array(5)).fill(1),
        invert: (new Array(5)).fill(false),
        position: (new Array(5 * 3)).fill(0),
        rotation: (new Array(5 * 4)).fill(0),
        scale: (new Array(5 * 3)).fill(1),
    };
    for (let i = 0, il = props.objects.length; i < il; ++i) {
        const p = props.objects[i];
        type[i] = Clipping.Type[p.type];
        invert[i] = p.invert;
        Vec3.toArray(p.position, position, i * 3);
        Quat.toArray(Quat.setAxisAngle(tmpQuat, p.rotation.axis, degToRad(p.rotation.angle)), rotation, i * 4);
        Vec3.toArray(p.scale, scale, i * 3);
    }
    return {
        variant: props.variant,
        objects: { count: props.objects.length, type, invert, position, rotation, scale }
    };
}

namespace Renderer {
    export function create(ctx: WebGLContext, props: Partial<RendererProps> = {}): Renderer {
        const { gl, state, stats, extensions: { fragDepth } } = ctx;
        const p = PD.merge(RendererParams, PD.getDefaultValues(RendererParams), props);
        const light = getLight(p.light);
        const clip = getClip(p.clip);

        const viewport = Viewport();
        const drawingBufferSize = Vec2.create(gl.drawingBufferWidth, gl.drawingBufferHeight);
        const bgColor = Color.toVec3Normalized(Vec3(), p.backgroundColor);

        let transparentBackground = false;

        const nullDepthTexture = createNullTexture(gl);
        const sharedTexturesList: Textures = [
            ['tDepth', nullDepthTexture]
        ];

        const view = Mat4();
        const invView = Mat4();
        const modelView = Mat4();
        const invModelView = Mat4();
        const invProjection = Mat4();
        const modelViewProjection = Mat4();
        const invModelViewProjection = Mat4();

        const cameraDir = Vec3();
        const viewOffset = Vec2();

        const ambientColor = Vec3();
        Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);

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
            uViewport: ValueCell.create(Viewport.toVec4(Vec4(), viewport)),
            uDrawingBufferSize: ValueCell.create(drawingBufferSize),

            uCameraPosition: ValueCell.create(Vec3()),
            uCameraDir: ValueCell.create(cameraDir),
            uNear: ValueCell.create(1),
            uFar: ValueCell.create(10000),
            uFogNear: ValueCell.create(1),
            uFogFar: ValueCell.create(10000),
            uFogColor: ValueCell.create(bgColor),

            uRenderWboit: ValueCell.create(false),
            uMarkingDepthTest: ValueCell.create(false),

            uTransparentBackground: ValueCell.create(false),

            uClipObjectType: ValueCell.create(clip.objects.type),
            uClipObjectInvert: ValueCell.create(clip.objects.invert),
            uClipObjectPosition: ValueCell.create(clip.objects.position),
            uClipObjectRotation: ValueCell.create(clip.objects.rotation),
            uClipObjectScale: ValueCell.create(clip.objects.scale),

            uLightDirection: ValueCell.create(light.direction),
            uLightColor: ValueCell.create(light.color),
            uAmbientColor: ValueCell.create(ambientColor),

            uPickingAlphaThreshold: ValueCell.create(p.pickingAlphaThreshold),

            uInteriorDarkening: ValueCell.create(p.interiorDarkening),
            uInteriorColorFlag: ValueCell.create(p.interiorColorFlag),
            uInteriorColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.interiorColor)),

            uHighlightColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.highlightColor)),
            uSelectColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.selectColor)),
            uHighlightStrength: ValueCell.create(p.highlightStrength),
            uSelectStrength: ValueCell.create(p.selectStrength),
            uMarkerPriority: ValueCell.create(p.markerPriority),

            uXrayEdgeFalloff: ValueCell.create(p.xrayEdgeFalloff),
        };
        const globalUniformList = Object.entries(globalUniforms);

        let globalUniformsNeedUpdate = true;

        const renderObject = (r: GraphicsRenderable, variant: GraphicsRenderVariant) => {
            if (r.state.disposed || !r.state.visible || (!r.state.pickable && variant[0] === 'p')) {
                return;
            }

            let definesNeedUpdate = false;
            if (r.state.noClip) {
                if (r.values.dClipObjectCount.ref.value !== 0) {
                    ValueCell.update(r.values.dClipObjectCount, 0);
                    definesNeedUpdate = true;
                }
            } else {
                if (r.values.dClipObjectCount.ref.value !== clip.objects.count) {
                    ValueCell.update(r.values.dClipObjectCount, clip.objects.count);
                    definesNeedUpdate = true;
                }
                if (r.values.dClipVariant.ref.value !== clip.variant) {
                    ValueCell.update(r.values.dClipVariant, clip.variant);
                    definesNeedUpdate = true;
                }
            }
            if (r.values.dLightCount.ref.value !== light.count) {
                ValueCell.update(r.values.dLightCount, light.count);
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

            if (r.values.dRenderMode) { // indicates direct-volume
                if ((variant[0] === 'p' || variant === 'depth') && r.values.dRenderMode.ref.value === 'volume') {
                    return; // no picking/depth in volume mode
                }

                // culling done in fragment shader
                state.disable(gl.CULL_FACE);
                state.frontFace(gl.CCW);

                if (variant === 'colorBlended') {
                    // depth test done manually in shader against `depthTexture`
                    // still need to enable when fragDepth can be used to write depth
                    if (r.values.dRenderMode.ref.value === 'volume' || !fragDepth) {
                        state.disable(gl.DEPTH_TEST);
                        state.depthMask(false);
                    } else {
                        state.enable(gl.DEPTH_TEST);
                        state.depthMask(r.values.uAlpha.ref.value === 1.0);
                    }
                }
            } else {
                if (r.values.dDoubleSided) {
                    if (r.values.dDoubleSided.ref.value || r.values.hasReflection.ref.value) {
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
            }

            r.render(variant, sharedTexturesList);
        };

        const update = (camera: ICamera) => {
            ValueCell.update(globalUniforms.uView, camera.view);
            ValueCell.update(globalUniforms.uInvView, Mat4.invert(invView, camera.view));
            ValueCell.update(globalUniforms.uProjection, camera.projection);
            ValueCell.update(globalUniforms.uInvProjection, Mat4.invert(invProjection, camera.projection));

            ValueCell.updateIfChanged(globalUniforms.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
            ValueCell.update(globalUniforms.uViewOffset, camera.viewOffset.enabled ? Vec2.set(viewOffset, camera.viewOffset.offsetX * 16, camera.viewOffset.offsetY * 16) : Vec2.set(viewOffset, 0, 0));

            ValueCell.update(globalUniforms.uCameraPosition, camera.state.position);
            ValueCell.update(globalUniforms.uCameraDir, Vec3.normalize(cameraDir, Vec3.sub(cameraDir, camera.state.target, camera.state.position)));

            ValueCell.updateIfChanged(globalUniforms.uFar, camera.far);
            ValueCell.updateIfChanged(globalUniforms.uNear, camera.near);
            ValueCell.updateIfChanged(globalUniforms.uFogFar, camera.fogFar);
            ValueCell.updateIfChanged(globalUniforms.uFogNear, camera.fogNear);
            ValueCell.updateIfChanged(globalUniforms.uTransparentBackground, transparentBackground);
        };

        const updateInternal = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null, renderWboit: boolean, markingDepthTest: boolean) => {
            arrayMapUpsert(sharedTexturesList, 'tDepth', depthTexture || nullDepthTexture);

            ValueCell.update(globalUniforms.uModel, group.view);
            ValueCell.update(globalUniforms.uModelView, Mat4.mul(modelView, group.view, camera.view));
            ValueCell.update(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView));
            ValueCell.update(globalUniforms.uModelViewProjection, Mat4.mul(modelViewProjection, modelView, camera.projection));
            ValueCell.update(globalUniforms.uInvModelViewProjection, Mat4.invert(invModelViewProjection, modelViewProjection));

            ValueCell.updateIfChanged(globalUniforms.uRenderWboit, renderWboit);
            ValueCell.updateIfChanged(globalUniforms.uMarkingDepthTest, markingDepthTest);

            state.enable(gl.SCISSOR_TEST);
            state.colorMask(true, true, true, true);

            const { x, y, width, height } = viewport;
            gl.viewport(x, y, width, height);
            gl.scissor(x, y, width, height);

            globalUniformsNeedUpdate = true;
            state.currentRenderItemId = -1;
        };

        const renderPick = (group: Scene.Group, camera: ICamera, variant: GraphicsRenderVariant, depthTexture: Texture | null) => {
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                if (!renderables[i].state.colorOnly) {
                    renderObject(renderables[i], variant);
                }
            }
        };

        const renderDepth = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                renderObject(renderables[i], 'depth');
            }
        };

        const renderMarkingDepth = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                if (r.values.markerAverage.ref.value !== 1) {
                    renderObject(renderables[i], 'markingDepth');
                }
            }
        };

        const renderMarkingMask = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, false, !!depthTexture);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                if (r.values.markerAverage.ref.value > 0) {
                    renderObject(renderables[i], 'markingMask');
                }
            }
        };

        const renderBlended = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            renderBlendedOpaque(group, camera, depthTexture);
            renderBlendedTransparent(group, camera, depthTexture);
        };

        const renderBlendedOpaque = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                if (r.state.opaque) {
                    renderObject(r, 'colorBlended');
                }
            }
        };

        const renderBlendedTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.enable(gl.DEPTH_TEST);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;

            if (transparentBackground) {
                state.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
            } else {
                state.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
            }
            state.enable(gl.BLEND);

            state.depthMask(true);
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                if (!r.state.opaque && r.state.writeDepth) {
                    renderObject(r, 'colorBlended');
                }
            }

            state.depthMask(false);
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                if (!r.state.opaque && !r.state.writeDepth) {
                    renderObject(r, 'colorBlended');
                }
            }
        };

        const renderBlendedVolumeOpaque = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
            state.enable(gl.BLEND);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                if (alpha === 1 && r.values.transparencyAverage.ref.value !== 1 && !r.values.dXrayShaded?.ref.value) {
                    renderObject(r, 'colorBlended');
                }
            }
        };

        const renderBlendedVolumeTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
            state.enable(gl.BLEND);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                if (alpha < 1 || r.values.transparencyAverage.ref.value > 0 || r.values.dXrayShaded?.ref.value) {
                    renderObject(r, 'colorBlended');
                }
            }
        };

        const renderWboitOpaque = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, false, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                if (alpha === 1 && r.values.transparencyAverage.ref.value !== 1 && r.values.dRenderMode?.ref.value !== 'volume' && r.values.dPointStyle?.ref.value !== 'fuzzy' && !r.values.dXrayShaded?.ref.value) {
                    renderObject(r, 'colorWboit');
                }
            }
        };

        const renderWboitTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            updateInternal(group, camera, depthTexture, true, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                if (alpha < 1 || r.values.transparencyAverage.ref.value > 0 || r.values.dRenderMode?.ref.value === 'volume' || r.values.dPointStyle?.ref.value === 'fuzzy' || !!r.values.uBackgroundColor || r.values.dXrayShaded?.ref.value) {
                    renderObject(r, 'colorWboit');
                }
            }
        };

        return {
            clear: (toBackgroundColor: boolean) => {
                state.enable(gl.SCISSOR_TEST);
                state.enable(gl.DEPTH_TEST);
                state.colorMask(true, true, true, true);
                state.depthMask(true);

                if (transparentBackground) {
                    state.clearColor(0, 0, 0, 0);
                } else if (toBackgroundColor) {
                    state.clearColor(bgColor[0], bgColor[1], bgColor[2], 1);
                } else {
                    state.clearColor(1, 1, 1, 1);
                }
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
            },
            clearDepth: () => {
                state.enable(gl.SCISSOR_TEST);
                state.enable(gl.DEPTH_TEST);
                state.depthMask(true);
                gl.clear(gl.DEPTH_BUFFER_BIT);
            },
            update,

            renderPick,
            renderDepth,
            renderMarkingDepth,
            renderMarkingMask,
            renderBlended,
            renderBlendedOpaque,
            renderBlendedTransparent,
            renderBlendedVolumeOpaque,
            renderBlendedVolumeTransparent,
            renderWboitOpaque,
            renderWboitTransparent,

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

                if (props.interiorDarkening !== undefined && props.interiorDarkening !== p.interiorDarkening) {
                    p.interiorDarkening = props.interiorDarkening;
                    ValueCell.update(globalUniforms.uInteriorDarkening, p.interiorDarkening);
                }
                if (props.interiorColorFlag !== undefined && props.interiorColorFlag !== p.interiorColorFlag) {
                    p.interiorColorFlag = props.interiorColorFlag;
                    ValueCell.update(globalUniforms.uInteriorColorFlag, p.interiorColorFlag);
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
                if (props.highlightStrength !== undefined && props.highlightStrength !== p.highlightStrength) {
                    p.highlightStrength = props.highlightStrength;
                    ValueCell.update(globalUniforms.uHighlightStrength, p.highlightStrength);
                }
                if (props.selectStrength !== undefined && props.selectStrength !== p.selectStrength) {
                    p.selectStrength = props.selectStrength;
                    ValueCell.update(globalUniforms.uSelectStrength, p.selectStrength);
                }
                if (props.markerPriority !== undefined && props.markerPriority !== p.markerPriority) {
                    p.markerPriority = props.markerPriority;
                    ValueCell.update(globalUniforms.uMarkerPriority, p.markerPriority);
                }

                if (props.xrayEdgeFalloff !== undefined && props.xrayEdgeFalloff !== p.xrayEdgeFalloff) {
                    p.xrayEdgeFalloff = props.xrayEdgeFalloff;
                    ValueCell.update(globalUniforms.uXrayEdgeFalloff, p.xrayEdgeFalloff);
                }

                if (props.light !== undefined && !deepEqual(props.light, p.light)) {
                    p.light = props.light;
                    Object.assign(light, getLight(props.light, light));
                    ValueCell.update(globalUniforms.uLightDirection, light.direction);
                    ValueCell.update(globalUniforms.uLightColor, light.color);
                }
                if (props.ambientColor !== undefined && props.ambientColor !== p.ambientColor) {
                    p.ambientColor = props.ambientColor;
                    Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);
                    ValueCell.update(globalUniforms.uAmbientColor, ambientColor);
                }
                if (props.ambientIntensity !== undefined && props.ambientIntensity !== p.ambientIntensity) {
                    p.ambientIntensity = props.ambientIntensity;
                    Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);
                    ValueCell.update(globalUniforms.uAmbientColor, ambientColor);
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
                gl.scissor(x, y, width, height);
                if (x !== viewport.x || y !== viewport.y || width !== viewport.width || height !== viewport.height) {
                    Viewport.set(viewport, x, y, width, height);
                    ValueCell.update(globalUniforms.uViewport, Vec4.set(globalUniforms.uViewport.ref.value, x, y, width, height));
                }
            },
            setTransparentBackground: (value: boolean) => {
                transparentBackground = value;
            },
            setDrawingBufferSize: (width: number, height: number) => {
                if (width !== drawingBufferSize[0] || height !== drawingBufferSize[1]) {
                    ValueCell.update(globalUniforms.uDrawingBufferSize, Vec2.set(drawingBufferSize, width, height));
                }
            },
            setPixelRatio: (value: number) => {
                ValueCell.update(globalUniforms.uPixelRatio, value);
            },

            props: p,
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

export { Renderer };