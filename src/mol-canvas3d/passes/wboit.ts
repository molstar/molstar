/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { evaluateWboit_frag } from '../../mol-gl/shader/evaluate-wboit.frag';
import { Framebuffer } from '../../mol-gl/webgl/framebuffer';
import { Vec2 } from '../../mol-math/linear-algebra';
import { isDebugMode, isTimingMode } from '../../mol-util/debug';
import { isWebGL2 } from '../../mol-gl/webgl/compat';
import { Renderbuffer } from '../../mol-gl/webgl/renderbuffer';

const EvaluateWboitSchema = {
    ...QuadSchema,
    tWboitA: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tWboitB: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    uTexSize: UniformSpec('v2'),
};
const EvaluateWboitShaderCode = ShaderCode('evaluate-wboit', quad_vert, evaluateWboit_frag);
type EvaluateWboitRenderable = ComputeRenderable<Values<typeof EvaluateWboitSchema>>

function getEvaluateWboitRenderable(ctx: WebGLContext, wboitATexture: Texture, wboitBTexture: Texture): EvaluateWboitRenderable {
    const values: Values<typeof EvaluateWboitSchema> = {
        ...QuadValues,
        tWboitA: ValueCell.create(wboitATexture),
        tWboitB: ValueCell.create(wboitBTexture),
        uTexSize: ValueCell.create(Vec2.create(wboitATexture.getWidth(), wboitATexture.getHeight())),
    };

    const schema = { ...EvaluateWboitSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', EvaluateWboitShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

export class WboitPass {
    private readonly renderable: EvaluateWboitRenderable;

    private readonly framebuffer: Framebuffer;
    private readonly textureA: Texture;
    private readonly textureB: Texture;
    private readonly depthRenderbuffer: Renderbuffer;

    private _supported = false;
    get supported() {
        return this._supported;
    }

    bind() {
        const { state, gl } = this.webgl;

        this.framebuffer.bind();

        state.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);

        state.disable(gl.DEPTH_TEST);

        state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ZERO, gl.ONE_MINUS_SRC_ALPHA);
        state.enable(gl.BLEND);
    }

    render() {
        if (isTimingMode) this.webgl.timer.mark('WboitPass.render');
        const { state, gl } = this.webgl;

        state.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
        state.enable(gl.BLEND);

        this.renderable.update();
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('WboitPass.render');
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;
        if (width !== w || height !== h) {
            this.textureA.define(width, height);
            this.textureB.define(width, height);
            this.depthRenderbuffer.setSize(width, height);
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
        }
    }

    reset() {
        if (this._supported) this._init();
    }

    private _init() {
        const { extensions: { drawBuffers } } = this.webgl;

        this.framebuffer.bind();
        drawBuffers!.drawBuffers([
            drawBuffers!.COLOR_ATTACHMENT0,
            drawBuffers!.COLOR_ATTACHMENT1,
        ]);

        this.textureA.attachFramebuffer(this.framebuffer, 'color0');
        this.textureB.attachFramebuffer(this.framebuffer, 'color1');

        this.depthRenderbuffer.attachFramebuffer(this.framebuffer);
    }

    static isSupported(webgl: WebGLContext) {
        const { extensions: { drawBuffers, textureFloat, colorBufferFloat, depthTexture } } = webgl;
        if (!textureFloat || !colorBufferFloat || !depthTexture || !drawBuffers) {
            if (isDebugMode) {
                const missing: string[] = [];
                if (!textureFloat) missing.push('textureFloat');
                if (!colorBufferFloat) missing.push('colorBufferFloat');
                if (!depthTexture) missing.push('depthTexture');
                if (!drawBuffers) missing.push('drawBuffers');
                console.log(`Missing "${missing.join('", "')}" extensions required for "wboit"`);
            }
            return false;
        } else {
            return true;
        }
    }

    constructor(private webgl: WebGLContext, width: number, height: number) {
        if (!WboitPass.isSupported(webgl)) return;

        const { resources, gl } = webgl;

        this.textureA = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        this.textureA.define(width, height);

        this.textureB = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        this.textureB.define(width, height);

        this.depthRenderbuffer = isWebGL2(gl)
            ? resources.renderbuffer('depth32f', 'depth', width, height)
            : resources.renderbuffer('depth16', 'depth', width, height);

        this.renderable = getEvaluateWboitRenderable(webgl, this.textureA, this.textureB);
        this.framebuffer = resources.framebuffer();

        this._supported = true;
        this._init();
    }
}