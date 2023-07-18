/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Viewport } from '../mol-canvas3d/camera/util';
import { ICamera } from '../mol-canvas3d/camera';
import { Scene } from './scene';
import { WebGLContext } from './webgl/context';
import { Mat4, Vec3, Vec4, Vec2 } from '../mol-math/linear-algebra';
import { GraphicsRenderable } from './renderable';
import { Color } from '../mol-util/color';
import { ValueCell, deepEqual } from '../mol-util';
import { GlobalUniformValues } from './renderable/schema';
import { GraphicsRenderVariant } from './webgl/render-item';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { degToRad } from '../mol-math/misc';
import { Texture, Textures } from './webgl/texture';
import { arrayMapUpsert } from '../mol-util/array';
import { clamp } from '../mol-math/interpolate';
import { isTimingMode } from '../mol-util/debug';

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

export enum PickType {
    None = 0,
    Object = 1,
    Instance = 2,
    Group = 3,
}

export enum MarkingType {
    None = 0,
    Depth = 1,
    Mask = 2,
}

interface Renderer {
    readonly stats: RendererStats
    readonly props: Readonly<RendererProps>
    readonly light: Readonly<Light>

    clear: (toBackgroundColor: boolean, ignoreTransparentBackground?: boolean) => void
    clearDepth: (packed?: boolean) => void
    update: (camera: ICamera, scene: Scene) => void

    renderPick: (group: Scene.Group, camera: ICamera, variant: 'pick' | 'depth', depthTexture: Texture | null, pickType: PickType) => void
    renderDepth: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderDepthOpaque: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderDepthTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderMarkingDepth: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderMarkingMask: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlended: (group: Scene, camera: ICamera) => void
    renderBlendedOpaque: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlendedTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderBlendedVolume: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderWboitOpaque: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderWboitTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderDpoitOpaque: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void
    renderDpoitTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null, dpoitTextures: { depth: Texture, frontColor: Texture, backColor: Texture }) => void
    renderDpoitVolume: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void

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

    colorMarker: PD.Boolean(true, { description: 'Enable color marker' }),
    highlightColor: PD.Color(Color.fromNormalizedRgb(1.0, 0.4, 0.6)),
    selectColor: PD.Color(Color.fromNormalizedRgb(0.2, 1.0, 0.1)),
    dimColor: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
    highlightStrength: PD.Numeric(0.3, { min: 0.0, max: 1.0, step: 0.1 }),
    selectStrength: PD.Numeric(0.3, { min: 0.0, max: 1.0, step: 0.1 }),
    dimStrength: PD.Numeric(0.0, { min: 0.0, max: 1.0, step: 0.1 }),
    markerPriority: PD.Select(1, [[1, 'Highlight'], [2, 'Select']]),

    xrayEdgeFalloff: PD.Numeric(1, { min: 0.0, max: 3.0, step: 0.1 }),
    exposure: PD.Numeric(1, { min: 0.0, max: 3.0, step: 0.01 }),

    light: PD.ObjectList({
        inclination: PD.Numeric(150, { min: 0, max: 180, step: 1 }),
        azimuth: PD.Numeric(320, { min: 0, max: 360, step: 1 }),
        color: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
        intensity: PD.Numeric(0.6, { min: 0.0, max: 5.0, step: 0.01 }),
    }, o => Color.toHexString(o.color), { defaultValue: [{
        inclination: 150,
        azimuth: 320,
        color: Color.fromNormalizedRgb(1.0, 1.0, 1.0),
        intensity: 0.6
    }] }),
    ambientColor: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
    ambientIntensity: PD.Numeric(0.4, { min: 0.0, max: 2.0, step: 0.01 }),
};
export type RendererProps = PD.Values<typeof RendererParams>

export type Light = {
    count: number
    direction: number[]
    color: number[]
}

const tmpDir = Vec3();
const tmpColor = Vec3();
function getLight(props: RendererProps['light'], light?: Light): Light {
    const count = props.length;
    const { direction, color } = light || {
        direction: (new Array(count * 3)).fill(0),
        color: (new Array(count * 3)).fill(0),
    };
    for (let i = 0; i < count; ++i) {
        const p = props[i];
        Vec3.directionFromSpherical(tmpDir, degToRad(p.inclination), degToRad(p.azimuth), 1);
        Vec3.toArray(tmpDir, direction, i * 3);
        Vec3.scale(tmpColor, Color.toVec3Normalized(tmpColor, p.color), p.intensity);
        Vec3.toArray(tmpColor, color, i * 3);
    }
    return { count, direction, color };
}

namespace Renderer {
    const enum Flag {
        None = 0,
        BlendedFront = 1,
        BlendedBack = 2,
    }

    const enum Mask {
        All = 0,
        Opaque = 1,
        Transparent = 2,
    }

    export function create(ctx: WebGLContext, props: Partial<RendererProps> = {}): Renderer {
        const { gl, state, stats } = ctx;
        const p = PD.merge(RendererParams, PD.getDefaultValues(RendererParams), props);
        const light = getLight(p.light);

        const viewport = Viewport();
        const drawingBufferSize = Vec2.create(gl.drawingBufferWidth, gl.drawingBufferHeight);
        const bgColor = Color.toVec3Normalized(Vec3(), p.backgroundColor);

        let transparentBackground = false;

        const emptyDepthTexture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        emptyDepthTexture.define(1, 1);
        emptyDepthTexture.load({ array: new Uint8Array([255, 255, 255, 255]), width: 1, height: 1 });
        const sharedTexturesList: Textures = [
            ['tDepth', emptyDepthTexture]
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
            uFog: ValueCell.create(true),
            uFogNear: ValueCell.create(1),
            uFogFar: ValueCell.create(10000),
            uFogColor: ValueCell.create(bgColor),

            uRenderMask: ValueCell.create(0),
            uMarkingDepthTest: ValueCell.create(false),
            uPickType: ValueCell.create(PickType.None),
            uMarkingType: ValueCell.create(MarkingType.None),

            uTransparentBackground: ValueCell.create(false),

            uLightDirection: ValueCell.create(light.direction),
            uLightColor: ValueCell.create(light.color),
            uAmbientColor: ValueCell.create(ambientColor),

            uPickingAlphaThreshold: ValueCell.create(p.pickingAlphaThreshold),

            uInteriorDarkening: ValueCell.create(p.interiorDarkening),
            uInteriorColorFlag: ValueCell.create(p.interiorColorFlag),
            uInteriorColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.interiorColor)),

            uHighlightColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.highlightColor)),
            uSelectColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.selectColor)),
            uDimColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.dimColor)),
            uHighlightStrength: ValueCell.create(p.highlightStrength),
            uSelectStrength: ValueCell.create(p.selectStrength),
            uDimStrength: ValueCell.create(p.dimStrength),
            uMarkerPriority: ValueCell.create(p.markerPriority),
            uMarkerAverage: ValueCell.create(0),

            uXrayEdgeFalloff: ValueCell.create(p.xrayEdgeFalloff),
            uExposure: ValueCell.create(p.exposure),
        };
        const globalUniformList = Object.entries(globalUniforms);

        let globalUniformsNeedUpdate = true;

        const renderObject = (r: GraphicsRenderable, variant: GraphicsRenderVariant, flag: Flag) => {
            if (r.state.disposed || !r.state.visible || (!r.state.pickable && variant === 'pick')) {
                return;
            }

            let definesNeedUpdate = false;
            if (r.values.dLightCount.ref.value !== light.count) {
                ValueCell.update(r.values.dLightCount, light.count);
                definesNeedUpdate = true;
            }
            if (r.values.dColorMarker.ref.value !== p.colorMarker) {
                ValueCell.update(r.values.dColorMarker, p.colorMarker);
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
                program.bindTextures(sharedTexturesList, 0);
                globalUniformsNeedUpdate = false;
            }

            if (r.values.dGeometryType.ref.value === 'directVolume') {
                if (variant !== 'colorDpoit' && variant !== 'colorWboit' && variant !== 'colorBlended') {
                    return; // only color supported
                }

                // culling done in fragment shader
                state.disable(gl.CULL_FACE);
                state.frontFace(gl.CCW);

                if (variant === 'colorBlended') {
                    // depth test done manually in shader against `depthTexture`
                    state.disable(gl.DEPTH_TEST);
                    state.depthMask(false);
                }
            } else if (flag === Flag.BlendedFront) {
                state.enable(gl.CULL_FACE);
                if (r.values.dFlipSided?.ref.value) {
                    state.frontFace(gl.CW);
                    state.cullFace(gl.FRONT);
                } else {
                    state.frontFace(gl.CCW);
                    state.cullFace(gl.BACK);
                }
            } else if (flag === Flag.BlendedBack) {
                state.enable(gl.CULL_FACE);
                if (r.values.dFlipSided?.ref.value) {
                    state.frontFace(gl.CW);
                    state.cullFace(gl.BACK);
                } else {
                    state.frontFace(gl.CCW);
                    state.cullFace(gl.FRONT);
                }
            } else {
                if (r.values.uDoubleSided) {
                    if (r.values.uDoubleSided.ref.value || r.values.hasReflection.ref.value) {
                        state.disable(gl.CULL_FACE);
                    } else {
                        state.enable(gl.CULL_FACE);
                    }
                } else {
                    // webgl default
                    state.disable(gl.CULL_FACE);
                }

                if (r.values.dFlipSided?.ref.value) {
                    state.frontFace(gl.CW);
                    state.cullFace(gl.FRONT);
                } else {
                    // webgl default
                    state.frontFace(gl.CCW);
                    state.cullFace(gl.BACK);
                }
            }

            r.render(variant, sharedTexturesList.length);
        };

        const update = (camera: ICamera, scene: Scene) => {
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

            ValueCell.updateIfChanged(globalUniforms.uMarkerAverage, scene.markerAverage);
        };

        const updateInternal = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null, renderMask: Mask, markingDepthTest: boolean) => {
            arrayMapUpsert(sharedTexturesList, 'tDepth', depthTexture || emptyDepthTexture);

            ValueCell.update(globalUniforms.uModel, group.view);
            ValueCell.update(globalUniforms.uModelView, Mat4.mul(modelView, camera.view, group.view));
            ValueCell.update(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView));
            ValueCell.update(globalUniforms.uModelViewProjection, Mat4.mul(modelViewProjection, modelView, camera.projection));
            ValueCell.update(globalUniforms.uInvModelViewProjection, Mat4.invert(invModelViewProjection, modelViewProjection));

            ValueCell.updateIfChanged(globalUniforms.uRenderMask, renderMask);
            ValueCell.updateIfChanged(globalUniforms.uMarkingDepthTest, markingDepthTest);

            state.enable(gl.SCISSOR_TEST);
            state.colorMask(true, true, true, true);

            const { x, y, width, height } = viewport;
            state.viewport(x, y, width, height);
            state.scissor(x, y, width, height);

            globalUniformsNeedUpdate = true;
            state.currentRenderItemId = -1;
        };

        const renderPick = (group: Scene.Group, camera: ICamera, variant: GraphicsRenderVariant, depthTexture: Texture | null, pickType: PickType) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderPick');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.All, false);
            ValueCell.updateIfChanged(globalUniforms.uPickType, pickType);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                if (!renderables[i].state.colorOnly) {
                    renderObject(renderables[i], variant, Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderPick');
        };

        const renderDepth = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderDepth');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.All, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                renderObject(renderables[i], 'depth', Flag.None);
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderDepth');
        };

        const renderDepthOpaque = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderDepthOpaque');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.Opaque, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
                if (r.state.opaque && r.values.transparencyAverage.ref.value !== 1 && !xrayShaded) {
                    renderObject(r, 'depth', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderDepthOpaque');
        };

        const renderDepthTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderDepthTransparent');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.Transparent, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
                if (!r.state.opaque || r.values.transparencyAverage.ref.value > 0 || xrayShaded) {
                    renderObject(r, 'depth', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderDepthTransparent');
        };

        const renderMarkingDepth = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderMarkingDepth');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.All, false);
            ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Depth);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                if (alpha !== 0 && r.values.markerAverage.ref.value !== 1) {
                    renderObject(renderables[i], 'marking', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderMarkingDepth');
        };

        const renderMarkingMask = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderMarkingMask');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.All, !!depthTexture);
            ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Mask);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                if (r.values.markerAverage.ref.value > 0) {
                    renderObject(renderables[i], 'marking', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderMarkingMask');
        };

        const renderBlended = (scene: Scene, camera: ICamera) => {
            if (scene.hasOpaque) {
                renderBlendedOpaque(scene, camera, null);
            }
            if (scene.opacityAverage < 1) {
                renderBlendedTransparent(scene, camera, null);
            }
        };

        const renderBlendedOpaque = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderBlendedOpaque');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.Opaque, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                if (r.state.opaque) {
                    renderObject(r, 'colorBlended', Flag.None);
                } else if (r.values.uDoubleSided?.ref.value && r.values.dTransparentBackfaces?.ref.value === 'opaque') {
                    renderObject(r, 'colorBlended', Flag.BlendedBack);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderBlendedOpaque');
        };

        const renderBlendedTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderBlendedTransparent');
            state.enable(gl.DEPTH_TEST);

            updateInternal(group, camera, depthTexture, Mask.Transparent, false);

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
                    renderObject(r, 'colorBlended', Flag.None);
                }
            }

            state.depthMask(false);
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                if ((!r.state.opaque && !r.state.writeDepth) || r.values.transparencyAverage.ref.value > 0) {
                    if (r.values.uDoubleSided?.ref.value) {
                        // render frontfaces and backfaces separately to avoid artefacts
                        if (r.values.dTransparentBackfaces?.ref.value !== 'opaque') {
                            renderObject(r, 'colorBlended', Flag.BlendedBack);
                        }
                        renderObject(r, 'colorBlended', Flag.BlendedFront);
                    } else {
                        renderObject(r, 'colorBlended', Flag.None);
                    }
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderBlendedTransparent');
        };

        const renderBlendedVolume = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderBlendedVolume');
            state.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
            state.enable(gl.BLEND);

            updateInternal(group, camera, depthTexture, Mask.Transparent, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                if (r.values.dGeometryType.ref.value === 'directVolume') {
                    renderObject(r, 'colorBlended', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderBlendedVolume');
        };

        const renderWboitOpaque = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderWboitOpaque');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.Opaque, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
                if ((alpha === 1 && r.values.transparencyAverage.ref.value !== 1 && r.values.dGeometryType.ref.value !== 'directVolume' && r.values.dPointStyle?.ref.value !== 'fuzzy' && !xrayShaded) || r.values.dTransparentBackfaces?.ref.value === 'opaque') {
                    renderObject(r, 'colorWboit', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderWboitOpaque');
        };

        const renderWboitTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderWboitTransparent');
            updateInternal(group, camera, depthTexture, Mask.Transparent, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
                if ((alpha < 1 && alpha !== 0) || r.values.transparencyAverage.ref.value > 0 || r.values.dGeometryType.ref.value === 'directVolume' || r.values.dPointStyle?.ref.value === 'fuzzy' || r.values.dGeometryType.ref.value === 'text' || xrayShaded) {
                    renderObject(r, 'colorWboit', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderWboitTransparent');
        };

        const renderDpoitOpaque = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderDpoitOpaque');
            state.disable(gl.BLEND);
            state.enable(gl.DEPTH_TEST);
            state.depthMask(true);

            updateInternal(group, camera, depthTexture, Mask.Opaque, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
                if ((alpha === 1 && r.values.transparencyAverage.ref.value !== 1 && r.values.dPointStyle?.ref.value !== 'fuzzy' && !xrayShaded) || r.values.dTransparentBackfaces?.ref.value === 'opaque') {
                    renderObject(r, 'colorDpoit', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderDpoitOpaque');
        };

        const renderDpoitTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null, dpoitTextures: { depth: Texture, frontColor: Texture, backColor: Texture }) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderDpoitTransparent');

            state.enable(gl.BLEND);

            arrayMapUpsert(sharedTexturesList, 'tDpoitDepth', dpoitTextures.depth);
            arrayMapUpsert(sharedTexturesList, 'tDpoitFrontColor', dpoitTextures.frontColor);
            arrayMapUpsert(sharedTexturesList, 'tDpoitBackColor', dpoitTextures.backColor);

            updateInternal(group, camera, depthTexture, Mask.Transparent, false);

            const { renderables } = group;

            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];

                // TODO: simplify, handle in renderable.state???
                // uAlpha is updated in "render" so we need to recompute it here
                const alpha = clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
                const xrayShaded = r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
                if ((alpha < 1 && alpha !== 0) || r.values.transparencyAverage.ref.value > 0 || r.values.dPointStyle?.ref.value === 'fuzzy' || r.values.dGeometryType.ref.value === 'text' || xrayShaded) {
                    renderObject(r, 'colorDpoit', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderDpoitTransparent');
        };

        const renderDpoitVolume = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
            if (isTimingMode) ctx.timer.mark('Renderer.renderDpoitVolume');
            state.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
            state.enable(gl.BLEND);

            updateInternal(group, camera, depthTexture, Mask.Transparent, false);

            const { renderables } = group;
            for (let i = 0, il = renderables.length; i < il; ++i) {
                const r = renderables[i];
                if (r.values.dGeometryType.ref.value === 'directVolume') {
                    renderObject(r, 'colorDpoit', Flag.None);
                }
            }
            if (isTimingMode) ctx.timer.markEnd('Renderer.renderDpoitVolume');
        };

        return {
            clear: (toBackgroundColor: boolean, ignoreTransparentBackground?: boolean) => {
                state.enable(gl.SCISSOR_TEST);
                state.enable(gl.DEPTH_TEST);
                state.colorMask(true, true, true, true);
                state.depthMask(true);

                if (transparentBackground && !ignoreTransparentBackground) {
                    state.clearColor(0, 0, 0, 0);
                } else if (toBackgroundColor) {
                    state.clearColor(bgColor[0], bgColor[1], bgColor[2], 1);
                } else {
                    state.clearColor(1, 1, 1, 1);
                }
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
            },
            clearDepth: (packed = false) => {
                state.enable(gl.SCISSOR_TEST);

                if (packed) {
                    state.colorMask(true, true, true, true);
                    state.clearColor(1, 1, 1, 1);
                    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
                } else {
                    state.enable(gl.DEPTH_TEST);
                    state.depthMask(true);
                    gl.clear(gl.DEPTH_BUFFER_BIT);
                }
            },
            update,

            renderPick,
            renderDepth,
            renderDepthOpaque,
            renderDepthTransparent,
            renderMarkingDepth,
            renderMarkingMask,
            renderBlended,
            renderBlendedOpaque,
            renderBlendedTransparent,
            renderBlendedVolume,
            renderWboitOpaque,
            renderWboitTransparent,
            renderDpoitOpaque,
            renderDpoitTransparent,
            renderDpoitVolume,

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

                if (props.colorMarker !== undefined && props.colorMarker !== p.colorMarker) {
                    p.colorMarker = props.colorMarker;
                }
                if (props.highlightColor !== undefined && props.highlightColor !== p.highlightColor) {
                    p.highlightColor = props.highlightColor;
                    ValueCell.update(globalUniforms.uHighlightColor, Color.toVec3Normalized(globalUniforms.uHighlightColor.ref.value, p.highlightColor));
                }
                if (props.selectColor !== undefined && props.selectColor !== p.selectColor) {
                    p.selectColor = props.selectColor;
                    ValueCell.update(globalUniforms.uSelectColor, Color.toVec3Normalized(globalUniforms.uSelectColor.ref.value, p.selectColor));
                }
                if (props.dimColor !== undefined && props.dimColor !== p.dimColor) {
                    p.dimColor = props.dimColor;
                    ValueCell.update(globalUniforms.uDimColor, Color.toVec3Normalized(globalUniforms.uDimColor.ref.value, p.dimColor));
                }
                if (props.highlightStrength !== undefined && props.highlightStrength !== p.highlightStrength) {
                    p.highlightStrength = props.highlightStrength;
                    ValueCell.update(globalUniforms.uHighlightStrength, p.highlightStrength);
                }
                if (props.selectStrength !== undefined && props.selectStrength !== p.selectStrength) {
                    p.selectStrength = props.selectStrength;
                    ValueCell.update(globalUniforms.uSelectStrength, p.selectStrength);
                }
                if (props.dimStrength !== undefined && props.dimStrength !== p.dimStrength) {
                    p.dimStrength = props.dimStrength;
                    ValueCell.update(globalUniforms.uDimStrength, p.dimStrength);
                }
                if (props.markerPriority !== undefined && props.markerPriority !== p.markerPriority) {
                    p.markerPriority = props.markerPriority;
                    ValueCell.update(globalUniforms.uMarkerPriority, p.markerPriority);
                }

                if (props.xrayEdgeFalloff !== undefined && props.xrayEdgeFalloff !== p.xrayEdgeFalloff) {
                    p.xrayEdgeFalloff = props.xrayEdgeFalloff;
                    ValueCell.update(globalUniforms.uXrayEdgeFalloff, p.xrayEdgeFalloff);
                }
                if (props.exposure !== undefined && props.exposure !== p.exposure) {
                    p.exposure = props.exposure;
                    ValueCell.update(globalUniforms.uExposure, p.exposure);
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
            },
            setViewport: (x: number, y: number, width: number, height: number) => {
                state.viewport(x, y, width, height);
                state.scissor(x, y, width, height);
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
            get light(): Light {
                return light;
            },
            dispose: () => {
                // TODO
            }
        };
    }
}

export { Renderer };
