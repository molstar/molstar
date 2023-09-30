/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { DefineSpec, TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture } from '../../mol-gl/webgl/texture';
import { Vec2 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { isTimingMode } from '../../mol-util/debug';
import { cas_frag } from '../../mol-gl/shader/cas.frag';

export const CasParams = {
    sharpness: PD.Numeric(0.5, { min: 0, max: 1, step: 0.05 }),
    denoise: PD.Boolean(true),
};
export type CasProps = PD.Values<typeof CasParams>

export class CasPass {
    private readonly renderable: CasRenderable;

    constructor(private webgl: WebGLContext, input: Texture) {
        this.renderable = getCasRenderable(webgl, input);
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
        ValueCell.update(this.renderable.values.uTexSizeInv, Vec2.set(this.renderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
    }

    update(input: Texture, props: CasProps) {
        const { values } = this.renderable;
        const { sharpness, denoise } = props;

        let needsUpdate = false;

        if (values.tColor.ref.value !== input) {
            ValueCell.update(this.renderable.values.tColor, input);
            needsUpdate = true;
        }

        ValueCell.updateIfChanged(values.uSharpness, 2 - 2 * Math.pow(sharpness, 0.25));

        if (values.dDenoise.ref.value !== denoise) needsUpdate = true;
        ValueCell.updateIfChanged(values.dDenoise, denoise);

        if (needsUpdate) {
            this.renderable.update();
        }
    }

    render(viewport: Viewport, target: RenderTarget | undefined) {
        if (isTimingMode) this.webgl.timer.mark('CasPass.render');
        if (target) {
            target.bind();
        } else {
            this.webgl.unbindFramebuffer();
        }
        this.updateState(viewport);
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('CasPass.render');
    }
}

//

const CasSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),

    uSharpness: UniformSpec('f'),
    dDenoise: DefineSpec('boolean'),
};
const CasShaderCode = ShaderCode('cas', quad_vert, cas_frag);
type CasRenderable = ComputeRenderable<Values<typeof CasSchema>>

function getCasRenderable(ctx: WebGLContext, colorTexture: Texture): CasRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof CasSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

        uSharpness: ValueCell.create(0.5),
        dDenoise: ValueCell.create(true),
    };

    const schema = { ...CasSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', CasShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}