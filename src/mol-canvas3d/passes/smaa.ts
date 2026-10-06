/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getAreaTexture, getSearchTexture } from './smaa-textures';
import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { loadImageTexture, Texture } from '../../mol-gl/webgl/texture';
import { Vec2, Vec4 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { blend_vert } from '../../mol-gl/shader/smaa/blend.vert';
import { blend_frag } from '../../mol-gl/shader/smaa/blend.frag';
import { weights_vert } from '../../mol-gl/shader/smaa/weights.vert';
import { weights_frag } from '../../mol-gl/shader/smaa/weights.frag';
import { edges_vert } from '../../mol-gl/shader/smaa/edges.vert';
import { edges_frag } from '../../mol-gl/shader/smaa/edges.frag';
import { Viewport } from '../camera/util';
import { isDebugMode, isTimingMode } from '../../mol-util/debug';

export const SmaaParams = {
    edgeThreshold: PD.Numeric(0.1, { min: 0.05, max: 0.15, step: 0.01 }),
    maxSearchSteps: PD.Numeric(16, { min: 0, max: 32, step: 1 }),
};
export type SmaaProps = PD.Values<typeof SmaaParams>

export class SmaaPass {
    private readonly edgesTarget: RenderTarget;
    private readonly weightsTarget: RenderTarget;

    private readonly edgesRenderable: EdgesRenderable;
    private readonly weightsRenderable: WeightsRenderable;
    private readonly blendRenderable: BlendRenderable;

    private _supported = false;
    get supported() {
        return this._supported;
    }

    constructor(private webgl: WebGLContext, input: Texture) {
        if (typeof HTMLImageElement === 'undefined') {
            if (isDebugMode) console.log(`Missing "HTMLImageElement" required for "SMAA"`);
            return;
        }

        const width = input.getWidth();
        const height = input.getHeight();

        this.edgesTarget = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');
        this.weightsTarget = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');

        this.edgesRenderable = getEdgesRenderable(webgl, input);
        this.weightsRenderable = getWeightsRenderable(webgl, this.edgesTarget.texture);
        this.blendRenderable = getBlendRenderable(webgl, input, this.weightsTarget.texture);

        this._supported = true;
    }

    getByteCount() {
        if (!this.supported) return 0;
        return this.edgesTarget.getByteCount() + this.weightsTarget.getByteCount();
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

        state.colorMask(true, true, true, true);
        state.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);

        ValueCell.update(this.edgesRenderable.values.uViewport, Viewport.toVec4(this.edgesRenderable.values.uViewport.ref.value, viewport));
        ValueCell.update(this.weightsRenderable.values.uViewport, Viewport.toVec4(this.weightsRenderable.values.uViewport.ref.value, viewport));
        ValueCell.update(this.blendRenderable.values.uViewport, Viewport.toVec4(this.blendRenderable.values.uViewport.ref.value, viewport));
    }

    setSize(width: number, height: number) {
        const w = this.edgesTarget.getWidth();
        const h = this.edgesTarget.getHeight();

        if (width !== w || height !== h) {
            this.edgesTarget.setSize(width, height);
            this.weightsTarget.setSize(width, height);

            ValueCell.update(this.edgesRenderable.values.uTexSizeInv, Vec2.set(this.edgesRenderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.weightsRenderable.values.uTexSizeInv, Vec2.set(this.weightsRenderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.blendRenderable.values.uTexSizeInv, Vec2.set(this.blendRenderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
        }
    }

    update(input: Texture, props: SmaaProps) {
        let edgesNeedsUpdate = false;
        if (this.edgesRenderable.values.tColor.ref.value !== input) {
            ValueCell.update(this.edgesRenderable.values.tColor, input);
            edgesNeedsUpdate = true;
        }
        if (this.edgesRenderable.values.dEdgeThreshold.ref.value !== props.edgeThreshold) {
            ValueCell.update(this.edgesRenderable.values.dEdgeThreshold, props.edgeThreshold);
            edgesNeedsUpdate = true;
        }
        if (edgesNeedsUpdate) this.edgesRenderable.update();

        if (this.weightsRenderable.values.dMaxSearchSteps.ref.value !== props.maxSearchSteps) {
            ValueCell.update(this.weightsRenderable.values.dMaxSearchSteps, props.maxSearchSteps);
            this.weightsRenderable.update();
        }

        if (this.blendRenderable.values.tColor.ref.value !== input) {
            ValueCell.update(this.blendRenderable.values.tColor, input);
            this.blendRenderable.update();
        }
    }

    render(viewport: Viewport, target: RenderTarget | undefined) {
        if (isTimingMode) this.webgl.timer.mark('SmaaPass.render');
        this.edgesTarget.bind();
        this.updateState(viewport);
        this.edgesRenderable.render();

        this.weightsTarget.bind();
        this.updateState(viewport);
        this.weightsRenderable.render();

        if (target) {
            target.bind();
        } else {
            this.webgl.bindDrawingBuffer();
        }
        this.updateState(viewport);
        this.blendRenderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('SmaaPass.render');
    }
}

//

const EdgesSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    uViewport: UniformSpec('v4'),

    dEdgeThreshold: DefineSpec('number')
};
const EdgesShaderCode = ShaderCode('smaa-edges', edges_vert, edges_frag);
type EdgesRenderable = ComputeRenderable<Values<typeof EdgesSchema>>

function getEdgesRenderable(ctx: WebGLContext, colorTexture: Texture): EdgesRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof EdgesSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uViewport: ValueCell.create(Vec4()),

        dEdgeThreshold: ValueCell.create(0.1),
    };

    const schema = { ...EdgesSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', EdgesShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

const WeightsSchema = {
    ...QuadSchema,
    tEdges: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tArea: TextureSpec('texture', 'rgb', 'ubyte', 'linear'),
    tSearch: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSizeInv: UniformSpec('v2'),
    uViewport: UniformSpec('v4'),

    dMaxSearchSteps: DefineSpec('number'),
};
const WeightsShaderCode = ShaderCode('smaa-weights', weights_vert, weights_frag);
type WeightsRenderable = ComputeRenderable<Values<typeof WeightsSchema>>

function getWeightsRenderable(ctx: WebGLContext, edgesTexture: Texture): WeightsRenderable {
    const width = edgesTexture.getWidth();
    const height = edgesTexture.getHeight();

    const areaTexture = ctx.resources.texture('image-uint8', 'rgb', 'ubyte', 'linear');
    const searchTexture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');

    const values: Values<typeof WeightsSchema> = {
        ...QuadValues,
        tEdges: ValueCell.create(edgesTexture),
        tArea: ValueCell.create(areaTexture),
        tSearch: ValueCell.create(searchTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uViewport: ValueCell.create(Vec4()),

        dMaxSearchSteps: ValueCell.create(16),
    };

    // Note: loading image textures requires `HTMLImageElement` to be available
    loadImageTexture(getAreaTexture(), values.tArea, areaTexture);
    loadImageTexture(getSearchTexture(), values.tSearch, searchTexture);

    const schema = { ...WeightsSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', WeightsShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

const BlendSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tWeights: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    uViewport: UniformSpec('v4'),
};
const BlendShaderCode = ShaderCode('smaa-blend', blend_vert, blend_frag);
type BlendRenderable = ComputeRenderable<Values<typeof BlendSchema>>

function getBlendRenderable(ctx: WebGLContext, colorTexture: Texture, weightsTexture: Texture): BlendRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof BlendSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        tWeights: ValueCell.create(weightsTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uViewport: ValueCell.create(Vec4()),
    };

    const schema = { ...BlendSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', BlendShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}
