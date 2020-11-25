/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from '../../mol-gl/renderable';
import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { DrawPass } from './draw';
import { Camera, ICamera } from '../../mol-canvas3d/camera';
import quad_vert from '../../mol-gl/shader/quad.vert';
import postprocessing_frag from '../../mol-gl/shader/postprocessing.frag';
import fxaa_frag from '../../mol-gl/shader/fxaa.frag';
import { StereoCamera } from '../camera/stereo';

const PostprocessingSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tPackedDepth: TextureSpec('texture', 'depth', 'ushort', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uFogNear: UniformSpec('f'),
    uFogFar: UniformSpec('f'),
    uFogColor: UniformSpec('v3'),

    dOcclusionEnable: DefineSpec('boolean'),
    dOcclusionKernelSize: DefineSpec('number'),
    uOcclusionBias: UniformSpec('f'),
    uOcclusionRadius: UniformSpec('f'),

    dOutlineEnable: DefineSpec('boolean'),
    uOutlineScale: UniformSpec('f'),
    uOutlineThreshold: UniformSpec('f'),
};
const PostprocessingShaderCode = ShaderCode('postprocessing', quad_vert, postprocessing_frag);
type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>

function getPostprocessingRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture): PostprocessingRenderable {
    const values: Values<typeof PostprocessingSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        tPackedDepth: ValueCell.create(depthTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uFogNear: ValueCell.create(10000),
        uFogFar: ValueCell.create(10000),
        uFogColor: ValueCell.create(Vec3.create(1, 1, 1)),

        dOcclusionEnable: ValueCell.create(false),
        dOcclusionKernelSize: ValueCell.create(4),
        uOcclusionBias: ValueCell.create(0.5),
        uOcclusionRadius: ValueCell.create(64),

        dOutlineEnable: ValueCell.create(false),
        uOutlineScale: ValueCell.create(1 * ctx.pixelRatio),
        uOutlineThreshold: ValueCell.create(0.8),
    };

    const schema = { ...PostprocessingSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', PostprocessingShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export const PostprocessingParams = {
    occlusion: PD.MappedStatic('off', {
        on: PD.Group({
            kernelSize: PD.Numeric(4, { min: 1, max: 32, step: 1 }),
            bias: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
            radius: PD.Numeric(64, { min: 0, max: 256, step: 1 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Darken occluded crevices with the ambient occlusion effect' }),
    outline: PD.MappedStatic('off', {
        on: PD.Group({
            scale: PD.Numeric(1, { min: 0, max: 10, step: 1 }),
            threshold: PD.Numeric(0.8, { min: 0, max: 5, step: 0.01 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Draw outline around 3D objects' }),
    antialiasing: PD.Boolean(true, { description: 'Fast Approximate Anti-Aliasing (FXAA)' })
};
export type PostprocessingProps = PD.Values<typeof PostprocessingParams>

export class PostprocessingPass {
    static isEnabled(props: PostprocessingProps) {
        return props.occlusion.name === 'on' || props.outline.name === 'on';
    }

    readonly target: RenderTarget

    private readonly tmpTarget: RenderTarget
    private readonly renderable: PostprocessingRenderable
    private readonly fxaa: FxaaRenderable

    constructor(private webgl: WebGLContext, private drawPass: DrawPass) {
        const { colorTarget, depthTexture } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.target = webgl.createRenderTarget(width, height, false);
        this.tmpTarget = webgl.createRenderTarget(width, height, false);
        this.renderable = getPostprocessingRenderable(webgl, colorTarget.texture, depthTexture);
        this.fxaa = getFxaaRenderable(webgl, this.tmpTarget.texture);
    }

    syncSize() {
        const width = this.drawPass.colorTarget.getWidth();
        const height = this.drawPass.colorTarget.getHeight();

        const [w, h] = this.renderable.values.uTexSize.ref.value;
        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            this.tmpTarget.setSize(width, height);
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.fxaa.values.uTexSize, Vec2.set(this.fxaa.values.uTexSize.ref.value, width, height));
        }
    }

    private updateState(camera: ICamera) {
        const { gl, state } = this.webgl;

        state.disable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = camera.viewport;
        gl.viewport(x, y, width, height);
        gl.scissor(x, y, width, height);

        state.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);
    }

    private _renderPostprocessing(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        const { values } = this.renderable;

        ValueCell.updateIfChanged(values.uFar, camera.far);
        ValueCell.updateIfChanged(values.uNear, camera.near);
        ValueCell.updateIfChanged(values.uFogFar, camera.fogFar);
        ValueCell.updateIfChanged(values.uFogNear, camera.fogNear);

        let needsUpdate = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        if (values.dOrthographic.ref.value !== orthographic) needsUpdate = true;
        ValueCell.updateIfChanged(values.dOrthographic, orthographic);

        const occlusion = props.occlusion.name === 'on';
        if (values.dOcclusionEnable.ref.value !== occlusion) needsUpdate = true;
        ValueCell.updateIfChanged(this.renderable.values.dOcclusionEnable, occlusion);
        if (props.occlusion.name === 'on') {
            const { kernelSize } = props.occlusion.params;
            if (values.dOcclusionKernelSize.ref.value !== kernelSize) needsUpdate = true;
            ValueCell.updateIfChanged(values.dOcclusionKernelSize, kernelSize);
            ValueCell.updateIfChanged(values.uOcclusionBias, props.occlusion.params.bias);
            ValueCell.updateIfChanged(values.uOcclusionRadius, props.occlusion.params.radius);
        }

        const outline = props.outline.name === 'on';
        if (values.dOutlineEnable.ref.value !== outline) needsUpdate = true;
        ValueCell.updateIfChanged(values.dOutlineEnable, outline);
        if (props.outline.name === 'on') {
            ValueCell.updateIfChanged(values.uOutlineScale, props.outline.params.scale * this.webgl.pixelRatio);
            ValueCell.updateIfChanged(values.uOutlineThreshold, props.outline.params.threshold);
        }

        if (needsUpdate) {
            this.renderable.update();
        }

        if (props.antialiasing) {
            this.tmpTarget.bind();
        } else if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
        } else {
            this.target.bind();
        }

        this.updateState(camera);
        this.renderable.render();
    }

    private _renderFxaa(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        const input = (props.occlusion.name === 'on' || props.outline.name === 'on')
            ? this.tmpTarget.texture : this.drawPass.colorTarget.texture;
        if (this.fxaa.values.tColor.ref.value !== input) {
            ValueCell.update(this.fxaa.values.tColor, input);
            this.fxaa.update();
        }

        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
        } else {
            this.target.bind();
        }

        this.updateState(camera);
        this.fxaa.render();
    }

    private _render(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        if (props.occlusion.name === 'on' || props.outline.name === 'on' || !props.antialiasing) {
            this._renderPostprocessing(camera, toDrawingBuffer, props);
        }

        if (props.antialiasing) {
            this._renderFxaa(camera, toDrawingBuffer, props);
        }
    }

    render(camera: Camera | StereoCamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        if (StereoCamera.is(camera)) {
            this._render(camera.left, toDrawingBuffer, props);
            this._render(camera.right, toDrawingBuffer, props);
        } else {
            this._render(camera, toDrawingBuffer, props);
        }
    }
}

//

const FxaaSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),
};
const FxaaShaderCode = ShaderCode('fxaa', quad_vert, fxaa_frag);
type FxaaRenderable = ComputeRenderable<Values<typeof FxaaSchema>>

function getFxaaRenderable(ctx: WebGLContext, colorTexture: Texture): FxaaRenderable {
    const values: Values<typeof FxaaSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),
    };

    const schema = { ...FxaaSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', FxaaShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}
