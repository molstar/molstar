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
import { Vec2 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { fxaa_frag } from '../../mol-gl/shader/fxaa.frag';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { isTimingMode } from '../../mol-util/debug';

export const FxaaParams = {
    edgeThresholdMin: PD.Numeric(0.0312, { min: 0.0312, max: 0.0833, step: 0.0001 }, { description: 'Trims the algorithm from processing darks.' }),
    edgeThresholdMax: PD.Numeric(0.063, { min: 0.063, max: 0.333, step: 0.001 }, { description: 'The minimum amount of local contrast required to apply algorithm.' }),
    iterations: PD.Numeric(12, { min: 0, max: 16, step: 1 }, { description: 'Number of edge exploration steps.' }),
    subpixelQuality: PD.Numeric(0.30, { min: 0.00, max: 1.00, step: 0.01 }, { description: 'Choose the amount of sub-pixel aliasing removal.' }),
};
export type FxaaProps = PD.Values<typeof FxaaParams>

export class FxaaPass {
    private readonly renderable: FxaaRenderable;

    constructor(private webgl: WebGLContext, input: Texture) {
        this.renderable = getFxaaRenderable(webgl, input);
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

    update(input: Texture, props: FxaaProps) {
        const { values } = this.renderable;
        const { edgeThresholdMin, edgeThresholdMax, iterations, subpixelQuality } = props;

        let needsUpdate = false;

        if (values.tColor.ref.value !== input) {
            ValueCell.update(this.renderable.values.tColor, input);
            needsUpdate = true;
        }

        if (values.dEdgeThresholdMin.ref.value !== edgeThresholdMin) needsUpdate = true;
        ValueCell.updateIfChanged(values.dEdgeThresholdMin, edgeThresholdMin);

        if (values.dEdgeThresholdMax.ref.value !== edgeThresholdMax) needsUpdate = true;
        ValueCell.updateIfChanged(values.dEdgeThresholdMax, edgeThresholdMax);

        if (values.dIterations.ref.value !== iterations) needsUpdate = true;
        ValueCell.updateIfChanged(values.dIterations, iterations);

        if (values.dSubpixelQuality.ref.value !== subpixelQuality) needsUpdate = true;
        ValueCell.updateIfChanged(values.dSubpixelQuality, subpixelQuality);

        if (needsUpdate) {
            this.renderable.update();
        }
    }

    render(viewport: Viewport, target: RenderTarget | undefined) {
        if (isTimingMode) this.webgl.timer.mark('FxaaPass.render');
        if (target) {
            target.bind();
        } else {
            this.webgl.unbindFramebuffer();
        }
        this.updateState(viewport);
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('FxaaPass.render');
    }
}

//

const FxaaSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),

    dEdgeThresholdMin: DefineSpec('number'),
    dEdgeThresholdMax: DefineSpec('number'),
    dIterations: DefineSpec('number'),
    dSubpixelQuality: DefineSpec('number'),
};
const FxaaShaderCode = ShaderCode('fxaa', quad_vert, fxaa_frag);
type FxaaRenderable = ComputeRenderable<Values<typeof FxaaSchema>>

function getFxaaRenderable(ctx: WebGLContext, colorTexture: Texture): FxaaRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof FxaaSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

        dEdgeThresholdMin: ValueCell.create(0.0312),
        dEdgeThresholdMax: ValueCell.create(0.125),
        dIterations: ValueCell.create(12),
        dSubpixelQuality: ValueCell.create(0.3),
    };

    const schema = { ...FxaaSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', FxaaShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}