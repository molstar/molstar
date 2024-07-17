/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from '../../mol-gl/renderable';
import { Mat4, Vec2 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { DrawPass } from './draw';
import { ICamera } from '../../mol-canvas3d/camera';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { outlines_frag } from '../../mol-gl/shader/outlines.frag';
import { isTimingMode } from '../../mol-util/debug';
import { Color } from '../../mol-util/color';
import { PostprocessingProps } from './postprocessing';

export const OutlineParams = {
    scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }),
    threshold: PD.Numeric(0.33, { min: 0.01, max: 1, step: 0.01 }),
    color: PD.Color(Color(0x000000)),
    includeTransparent: PD.Boolean(true, { description: 'Whether to show outline for transparent objects' }),
};

export type OutlineProps = PD.Values<typeof OutlineParams>

export class OutlinePass {
    static isEnabled(props: PostprocessingProps) {
        return props.outline.name !== 'off';
    }

    readonly outlinesTarget: RenderTarget;
    private readonly outlinesRenderable: OutlinesRenderable;

    constructor(private readonly webgl: WebGLContext, readonly drawPass: DrawPass, width: number, height: number) {
        const { depthTextureTransparent, depthTextureOpaque } = drawPass;

        this.outlinesTarget = webgl.createRenderTarget(width, height, false);
        this.outlinesRenderable = getOutlinesRenderable(webgl, depthTextureOpaque, depthTextureTransparent, true);
    }

    setSize(width: number, height: number) {
        const [w, h] = this.outlinesRenderable.values.uTexSize.ref.value;
        if (width !== w || height !== h) {
            this.outlinesTarget.setSize(width, height);
            ValueCell.update(this.outlinesRenderable.values.uTexSize, Vec2.set(this.outlinesRenderable.values.uTexSize.ref.value, width, height));
        }
    }

    update(camera: ICamera, props: OutlineProps) {
        let needsUpdateOutlines = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        const transparentOutline = props.includeTransparent ?? true;
        const outlineThreshold = 50 * props.threshold * this.webgl.pixelRatio;

        ValueCell.updateIfChanged(this.outlinesRenderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.outlinesRenderable.values.uFar, camera.far);
        ValueCell.update(this.outlinesRenderable.values.uInvProjection, invProjection);
        if (this.outlinesRenderable.values.dTransparentOutline.ref.value !== transparentOutline) {
            needsUpdateOutlines = true;
            ValueCell.update(this.outlinesRenderable.values.dTransparentOutline, transparentOutline);
        }
        if (this.outlinesRenderable.values.dOrthographic.ref.value !== orthographic) {
            needsUpdateOutlines = true;
            ValueCell.update(this.outlinesRenderable.values.dOrthographic, orthographic);
        }
        ValueCell.updateIfChanged(this.outlinesRenderable.values.uOutlineThreshold, outlineThreshold);

        if (needsUpdateOutlines) {
            this.outlinesRenderable.update();
        }
    }

    render() {
        if (isTimingMode) this.webgl.timer.mark('OUTLINE.render');
        this.outlinesTarget.bind();
        this.outlinesRenderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('OUTLINE.render');
    }
}

export const OutlinesSchema = {
    ...QuadSchema,
    tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uInvProjection: UniformSpec('m4'),

    uOutlineThreshold: UniformSpec('f'),
    dTransparentOutline: DefineSpec('boolean'),
};
export type OutlinesRenderable = ComputeRenderable<Values<typeof OutlinesSchema>>

export function getOutlinesRenderable(ctx: WebGLContext, depthTextureOpaque: Texture, depthTextureTransparent: Texture, transparentOutline: boolean): OutlinesRenderable {
    const width = depthTextureOpaque.getWidth();
    const height = depthTextureOpaque.getHeight();

    const values: Values<typeof OutlinesSchema> = {
        ...QuadValues,
        tDepthOpaque: ValueCell.create(depthTextureOpaque),
        tDepthTransparent: ValueCell.create(depthTextureTransparent),
        uTexSize: ValueCell.create(Vec2.create(width, height)),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uInvProjection: ValueCell.create(Mat4.identity()),

        uOutlineThreshold: ValueCell.create(0.33),
        dTransparentOutline: ValueCell.create(transparentOutline),
    };

    const schema = { ...OutlinesSchema };
    const shaderCode = ShaderCode('outlines', quad_vert, outlines_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}