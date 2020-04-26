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
import { Camera } from '../../mol-canvas3d/camera';
import { produce } from 'immer';

import quad_vert from '../../mol-gl/shader/quad.vert';
import postprocessing_frag from '../../mol-gl/shader/postprocessing.frag';

const PostprocessingSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
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

    dPackedDepth: DefineSpec('boolean'),
};

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
    }, { cycle: true, description: 'Draw outline around 3D objects' })
};
export type PostprocessingProps = PD.Values<typeof PostprocessingParams>

type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>

function getPostprocessingRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture, packedDepth: boolean, props: Partial<PostprocessingProps>): PostprocessingRenderable {
    const p = { ...PD.getDefaultValues(PostprocessingParams), ...props };
    const values: Values<typeof PostprocessingSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        tDepth: ValueCell.create(depthTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uFogNear: ValueCell.create(10000),
        uFogFar: ValueCell.create(10000),
        uFogColor: ValueCell.create(Vec3.create(1, 1, 1)),

        dOcclusionEnable: ValueCell.create(p.occlusion.name === 'on'),
        dOcclusionKernelSize: ValueCell.create(p.occlusion.name === 'on' ? p.occlusion.params.kernelSize : 4),
        uOcclusionBias: ValueCell.create(p.occlusion.name === 'on' ? p.occlusion.params.bias : 0.5),
        uOcclusionRadius: ValueCell.create(p.occlusion.name === 'on' ? p.occlusion.params.radius : 64),

        dOutlineEnable: ValueCell.create(p.outline.name === 'on'),
        uOutlineScale: ValueCell.create((p.outline.name === 'on' ? p.outline.params.scale : 1) * ctx.pixelRatio),
        uOutlineThreshold: ValueCell.create(p.outline.name === 'on' ? p.outline.params.threshold : 0.8),

        dPackedDepth: ValueCell.create(packedDepth),
    };

    const schema = { ...PostprocessingSchema };
    const shaderCode = ShaderCode('postprocessing', quad_vert, postprocessing_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export class PostprocessingPass {
    target: RenderTarget
    props: PostprocessingProps
    renderable: PostprocessingRenderable

    constructor(private webgl: WebGLContext, private camera: Camera, drawPass: DrawPass, props: Partial<PostprocessingProps>) {
        const { gl } = webgl;
        this.target = webgl.createRenderTarget(gl.drawingBufferWidth, gl.drawingBufferHeight);
        this.props = { ...PD.getDefaultValues(PostprocessingParams), ...props };
        const { colorTarget, depthTexture, packedDepth } = drawPass;
        this.renderable = getPostprocessingRenderable(webgl, colorTarget.texture, depthTexture, packedDepth, this.props);
    }

    get enabled() {
        return this.props.occlusion.name === 'on' || this.props.outline.name === 'on';
    }

    setSize(width: number, height: number) {
        this.target.setSize(width, height);
        ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
    }

    setProps(props: Partial<PostprocessingProps>) {
        this.props = produce(this.props, p => {
            if (props.occlusion !== undefined) {
                p.occlusion.name = props.occlusion.name;
                ValueCell.updateIfChanged(this.renderable.values.dOcclusionEnable, props.occlusion.name === 'on');
                if (props.occlusion.name === 'on') {
                    p.occlusion.params = { ...props.occlusion.params };
                    ValueCell.updateIfChanged(this.renderable.values.dOcclusionKernelSize, props.occlusion.params.kernelSize);
                    ValueCell.updateIfChanged(this.renderable.values.uOcclusionBias, props.occlusion.params.bias);
                    ValueCell.updateIfChanged(this.renderable.values.uOcclusionRadius, props.occlusion.params.radius);
                }
            }

            if (props.outline !== undefined) {
                p.outline.name = props.outline.name;
                ValueCell.updateIfChanged(this.renderable.values.dOutlineEnable, props.outline.name === 'on');
                if (props.outline.name === 'on') {
                    p.outline.params = { ...props.outline.params };
                    ValueCell.updateIfChanged(this.renderable.values.uOutlineScale, props.outline.params.scale);
                    ValueCell.updateIfChanged(this.renderable.values.uOutlineThreshold, props.outline.params.threshold);
                }
            }
        });

        this.renderable.update();
    }

    render(toDrawingBuffer: boolean) {
        ValueCell.update(this.renderable.values.uFar, this.camera.far);
        ValueCell.update(this.renderable.values.uNear, this.camera.near);
        ValueCell.update(this.renderable.values.uFogFar, this.camera.fogFar);
        ValueCell.update(this.renderable.values.uFogNear, this.camera.fogNear);
        ValueCell.update(this.renderable.values.dOrthographic, this.camera.state.mode === 'orthographic' ? 1 : 0);

        const { gl, state } = this.webgl;
        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
            gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
        } else {
            this.target.bind();
        }
        state.disable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);
        this.renderable.render();
    }
}