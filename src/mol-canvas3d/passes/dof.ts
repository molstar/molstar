/**
 * Copyright (c) 2024-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture, createNullTexture } from '../../mol-gl/webgl/texture';
import { Mat4, Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { dof_frag } from '../../mol-gl/shader/dof.frag';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { isTimingMode } from '../../mol-util/debug';
import { ICamera } from '../../mol-canvas3d/camera';
import { Sphere3D } from '../../mol-math/geometry';
import { PostprocessingProps } from './postprocessing';

export const DofParams = {
    blurSize: PD.Numeric(9, { min: 1, max: 32, step: 1 }),
    blurSpread: PD.Numeric(1.0, { min: 0.0, max: 10.0, step: 0.1 }),
    inFocus: PD.Numeric(0.0, { min: -5000.0, max: 5000.0, step: 1.0 }, { description: 'Distance from the scene center that will be in focus' }),
    PPM: PD.Numeric(20.0, { min: 0.0, max: 5000.0, step: 0.1 }, { description: 'Size of the area that will be in focus' }),
    center: PD.Select('camera-target', PD.arrayToOptions(['scene-center', 'camera-target'])),
    mode: PD.Select('plane', PD.arrayToOptions(['plane', 'sphere'])),
};

export type DofProps = PD.Values<typeof DofParams>

export class DofPass {
    static isEnabled(props: PostprocessingProps) {
        return props.enabled && props.dof.name !== 'off';
    }

    readonly target: RenderTarget;
    private readonly renderable: DofRenderable;

    constructor(private webgl: WebGLContext, width: number, height: number) {
        this.target = webgl.createRenderTarget(width, height, false);

        const nullTexture = createNullTexture();
        this.renderable = getDofRenderable(webgl, nullTexture, nullTexture, nullTexture);
    }

    private updateState(viewport: Viewport) {
        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = viewport;
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        state.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);
    }

    setSize(width: number, height: number) {
        const w = this.target.texture.getWidth();
        const h = this.target.texture.getHeight();

        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
        }
    }

    update(camera: ICamera, input: Texture, depthOpaque: Texture, depthTransparent: Texture, props: DofProps, sphere: Sphere3D) {
        let needsUpdate = false;
        if (this.renderable.values.tColor.ref.value !== input) {
            ValueCell.update(this.renderable.values.tColor, input);
            needsUpdate = true;
        }
        if (this.renderable.values.tDepthOpaque.ref.value !== depthOpaque) {
            ValueCell.update(this.renderable.values.tDepthOpaque, depthOpaque);
            needsUpdate = true;
        }
        if (this.renderable.values.tDepthTransparent.ref.value !== depthTransparent) {
            ValueCell.update(this.renderable.values.tDepthTransparent, depthTransparent);
            needsUpdate = true;
        }
        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        const invProjection = this.renderable.values.uInvProjection.ref.value;
        Mat4.invert(invProjection, camera.projection);

        const [w, h] = this.renderable.values.uTexSize.ref.value;
        const v = camera.viewport;

        ValueCell.update(this.renderable.values.uProjection, camera.projection);
        ValueCell.update(this.renderable.values.uInvProjection, invProjection);
        ValueCell.update(this.renderable.values.uMode, props.mode === 'sphere' ? 1 : 0);

        Vec4.set(this.renderable.values.uBounds.ref.value,
            v.x / w,
            v.y / h,
            (v.x + v.width) / w,
            (v.y + v.height) / h
        );
        ValueCell.update(this.renderable.values.uBounds, this.renderable.values.uBounds.ref.value);

        ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
        ValueCell.updateIfChanged(this.renderable.values.dOrthographic, orthographic);

        const blurSize = Math.round(props.blurSize * this.webgl.pixelRatio);
        if (this.renderable.values.dBlurSize.ref.value !== blurSize) {
            ValueCell.update(this.renderable.values.dBlurSize, blurSize);
            needsUpdate = true;
        }

        const worldCenter = (props.center === 'scene-center' ? sphere.center : camera.state.target);
        const distance = Vec3.distance(camera.state.position, worldCenter);
        const inFocus = distance + props.inFocus;
        ValueCell.updateIfChanged(this.renderable.values.uInFocus, inFocus * camera.state.scale);

        // transform center in view space
        const center = this.renderable.values.uCenter.ref.value;
        Vec3.transformMat4(center, worldCenter, camera.view);
        ValueCell.update(this.renderable.values.uCenter, center);

        ValueCell.updateIfChanged(this.renderable.values.uBlurSpread, props.blurSpread);
        ValueCell.updateIfChanged(this.renderable.values.uPPM, props.PPM * camera.state.scale);

        if (needsUpdate) {
            this.renderable.update();
        }
    }

    render(viewport: Viewport, target: undefined | RenderTarget) {
        if (isTimingMode) this.webgl.timer.mark('DofPass.render');
        if (target) {
            target.bind();
        } else {
            this.webgl.bindDrawingBuffer();
        }
        this.updateState(viewport);
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('DofPass.render');
    }
}

//

const DofSchema = {
    ...QuadSchema,
    tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),
    uBounds: UniformSpec('v4'),
    uCenter: UniformSpec('v3'),
    uMode: UniformSpec('i'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),

    dBlurSize: DefineSpec('number'),
    uBlurSpread: UniformSpec('f'),
    uInFocus: UniformSpec('f'),
    uPPM: UniformSpec('f'),
};

const DofShaderCode = ShaderCode('dof', quad_vert, dof_frag);
type DofRenderable = ComputeRenderable<Values<typeof DofSchema>>

function getDofRenderable(ctx: WebGLContext, colorTexture: Texture, depthTextureOpaque: Texture, depthTextureTransparent: Texture): DofRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof DofSchema> = {
        ...QuadValues,
        tDepthOpaque: ValueCell.create(depthTextureOpaque),
        tDepthTransparent: ValueCell.create(depthTextureTransparent),
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(width, height)),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),
        uBounds: ValueCell.create(Vec4()),
        uCenter: ValueCell.create(Vec3()),
        uMode: ValueCell.create(0),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),

        dBlurSize: ValueCell.create(5),
        uBlurSpread: ValueCell.create(300.0),
        uInFocus: ValueCell.create(20.0),
        uPPM: ValueCell.create(20.0),
    };

    const schema = { ...DofSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', DofShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}