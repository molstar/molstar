/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture } from '../../mol-gl/webgl/texture';
import { Mat4, Vec2, Vec4 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { dof_frag } from '../../mol-gl/shader/dof.frag';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { isTimingMode } from '../../mol-util/debug';
import { ICamera } from '../../mol-canvas3d/camera';

export const DofParams = {
    blurSize: PD.Numeric(5, { min: 1, max: 32, step: 1 }),
    blurSpread: PD.Numeric(1.0, { min: 0.0, max: 10.0, step: 0.1 }),
    inFocus: PD.Numeric(100.0, { min: 0.0, max: 10000.0, step: 0.1 }),
    PPM: PD.Numeric(20.0, { min: 0.0, max: 5000.0, step: 0.1 }),
};

export type DofProps = PD.Values<typeof DofParams>

export class DofPass {
    private readonly renderable: DofRenderable;

    constructor(private webgl: WebGLContext, input: Texture, depth: Texture) {
        this.renderable = getDofRenderable(webgl, input, depth);
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
        ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
    }

    update(camera: ICamera, input: Texture, depth: Texture, props: DofProps) {
        let needsUpdate = false;
        if (this.renderable.values.tColor.ref.value !== input) {
            ValueCell.update(this.renderable.values.tColor, input);
            needsUpdate = true;
        }
        if (this.renderable.values.tDepth.ref.value !== depth) {
            ValueCell.update(this.renderable.values.tDepth, depth);
            needsUpdate = true;
        }
        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        const [w, h] = this.renderable.values.uTexSize.ref.value;
        const v = camera.viewport;

        ValueCell.update(this.renderable.values.uProjection, camera.projection);
        ValueCell.update(this.renderable.values.uInvProjection, invProjection);

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

        if (this.renderable.values.blurSize.ref.value !== props.blurSize) needsUpdate = true;
        ValueCell.update(this.renderable.values.blurSize, props.blurSize);
        if (this.renderable.values.blurSpread.ref.value !== props.blurSpread) needsUpdate = true;
        ValueCell.update(this.renderable.values.blurSpread, props.blurSpread);
        if (this.renderable.values.inFocus.ref.value !== props.inFocus) needsUpdate = true;
        ValueCell.update(this.renderable.values.inFocus, props.inFocus);
        if (this.renderable.values.PPM.ref.value !== props.PPM) needsUpdate = true;
        ValueCell.update(this.renderable.values.PPM, props.PPM);

        this.renderable.update();
        ValueCell.update(this.renderable.values.PPM, props.PPM);
        ValueCell.update(this.renderable.values.inFocus, props.inFocus);
        ValueCell.update(this.renderable.values.blurSpread, props.blurSpread);
        ValueCell.update(this.renderable.values.blurSize, props.blurSize);

        if (needsUpdate) {
            this.renderable.update();
        }
    }

    render(viewport: Viewport, target: RenderTarget | undefined) {
        if (isTimingMode) this.webgl.timer.mark('DofPass.render');
        if (target) {
            target.bind();
        } else {
            this.webgl.unbindFramebuffer();
        }
        this.updateState(viewport);
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('DofPass.render');
    }
}

//

const DofSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),
    uBounds: UniformSpec('v4'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),

    inFocus: UniformSpec('f'),
    PPM: UniformSpec('f'),
    blurSize: UniformSpec('i'),
    blurSpread: UniformSpec('f'),
};

const DofShaderCode = ShaderCode('dof', quad_vert, dof_frag);
type DofRenderable = ComputeRenderable<Values<typeof DofSchema>>

function getDofRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture): DofRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof DofSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(width, height)),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),
        uBounds: ValueCell.create(Vec4()),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),

        blurSize: ValueCell.create(5),
        blurSpread: ValueCell.create(300.0),
        inFocus: ValueCell.create(20.0),
        PPM: ValueCell.create(20.0),
    };

    const schema = { ...DofSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', DofShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}