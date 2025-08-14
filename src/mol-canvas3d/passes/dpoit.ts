/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from https://github.com/tsherif/webgl2examples, The MIT License, Copyright © 2017 Tarek Sherif, Shuai Shao
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
import { evaluateDpoit_frag } from '../../mol-gl/shader/evaluate-dpoit.frag';
import { blendBackDpoit_frag } from '../../mol-gl/shader/blend-back-dpoit.frag';
import { Framebuffer } from '../../mol-gl/webgl/framebuffer';
import { Vec2 } from '../../mol-math/linear-algebra';
import { isDebugMode, isTimingMode } from '../../mol-util/debug';
import { isWebGL2 } from '../../mol-gl/webgl/compat';

const BlendBackDpoitSchema = {
    ...QuadSchema,
    tDpoitBackColor: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    uTexSize: UniformSpec('v2'),
};
const BlendBackDpoitShaderCode = ShaderCode('blend-back-dpoit', quad_vert, blendBackDpoit_frag);
type BlendBackDpoitRenderable = ComputeRenderable<Values<typeof BlendBackDpoitSchema>>

function getBlendBackDpoitRenderable(ctx: WebGLContext, dopitBlendBackTexture: Texture): BlendBackDpoitRenderable {
    const values: Values<typeof BlendBackDpoitSchema> = {
        ...QuadValues,
        tDpoitBackColor: ValueCell.create(dopitBlendBackTexture),
        uTexSize: ValueCell.create(Vec2.create(dopitBlendBackTexture.getWidth(), dopitBlendBackTexture.getHeight())),
    };

    const schema = { ...BlendBackDpoitSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', BlendBackDpoitShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

const EvaluateDpoitSchema = {
    ...QuadSchema,
    tDpoitFrontColor: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    uTexSize: UniformSpec('v2'),
};
const EvaluateDpoitShaderCode = ShaderCode('evaluate-dpoit', quad_vert, evaluateDpoit_frag);
type EvaluateDpoitRenderable = ComputeRenderable<Values<typeof EvaluateDpoitSchema>>

function getEvaluateDpoitRenderable(ctx: WebGLContext, dpoitFrontColorTexture: Texture): EvaluateDpoitRenderable {
    const values: Values<typeof EvaluateDpoitSchema> = {
        ...QuadValues,
        tDpoitFrontColor: ValueCell.create(dpoitFrontColorTexture),
        uTexSize: ValueCell.create(Vec2.create(dpoitFrontColorTexture.getWidth(), dpoitFrontColorTexture.getHeight())),
    };

    const schema = { ...EvaluateDpoitSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', EvaluateDpoitShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export class DpoitPass {
    private readonly DEPTH_CLEAR_VALUE = -99999.0; // NOTE same constant is set in shaders
    private readonly MAX_DEPTH = 1.0;
    private readonly MIN_DEPTH = 0.0;

    private passCount = 0;
    private writeId: number;
    private readId: number;

    private readonly blendBackRenderable: BlendBackDpoitRenderable;
    private readonly renderable: EvaluateDpoitRenderable;

    private readonly depthFramebuffers: Framebuffer[];
    private readonly colorFramebuffers: Framebuffer[];

    private readonly depthTextures: Texture[];
    private readonly colorFrontTextures: Texture[];
    private readonly colorBackTextures: Texture[];

    private _supported = false;
    get supported() {
        return this._supported;
    }

    bind() {
        const { state, gl, extensions: { blendMinMax } } = this.webgl;

        // initialize
        this.passCount = 0;

        this.depthFramebuffers[0].bind();
        state.clearColor(this.DEPTH_CLEAR_VALUE, this.DEPTH_CLEAR_VALUE, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        this.depthFramebuffers[1].bind();
        state.clearColor(-this.MIN_DEPTH, this.MAX_DEPTH, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        this.colorFramebuffers[0].bind();
        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        this.colorFramebuffers[1].bind();
        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        this.depthFramebuffers[0].bind();
        state.blendEquation(blendMinMax!.MAX);
        state.depthMask(false);

        return {
            depth: this.depthTextures[1],
            frontColor: this.colorFrontTextures[1],
            backColor: this.colorBackTextures[1]
        };
    }

    bindDualDepthPeeling() {
        const { state, gl, extensions: { blendMinMax } } = this.webgl;

        this.readId = this.passCount % 2;
        this.writeId = 1 - this.readId; // ping-pong: 0 or 1

        this.passCount += 1; // increment for next pass

        this.depthFramebuffers[this.writeId].bind();
        state.clearColor(this.DEPTH_CLEAR_VALUE, this.DEPTH_CLEAR_VALUE, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        this.colorFramebuffers[this.writeId].bind();
        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        this.depthFramebuffers[this.writeId].bind();
        state.blendEquation(blendMinMax!.MAX);
        state.depthMask(false);

        return {
            depth: this.depthTextures[this.readId],
            frontColor: this.colorFrontTextures[this.readId],
            backColor: this.colorBackTextures[this.readId]
        };
    }

    renderBlendBack() {
        if (isTimingMode) this.webgl.timer.mark('DpoitPass.renderBlendBack');
        const { state, gl } = this.webgl;

        state.blendEquation(gl.FUNC_ADD);
        state.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);

        ValueCell.update(this.blendBackRenderable.values.tDpoitBackColor, this.colorBackTextures[this.writeId]);

        this.blendBackRenderable.update();
        this.blendBackRenderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('DpoitPass.renderBlendBack');
    }

    render() {
        if (isTimingMode) this.webgl.timer.mark('DpoitPass.render');
        const { state, gl } = this.webgl;

        state.blendEquation(gl.FUNC_ADD);
        state.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);

        ValueCell.update(this.renderable.values.tDpoitFrontColor, this.colorFrontTextures[this.writeId]);

        this.renderable.update();
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('DpoitPass.render');
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;
        if (width !== w || height !== h) {
            for (let i = 0; i < 2; i++) {
                this.depthTextures[i].define(width, height);
                this.colorFrontTextures[i].define(width, height);
                this.colorBackTextures[i].define(width, height);
            }
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.blendBackRenderable.values.uTexSize, Vec2.set(this.blendBackRenderable.values.uTexSize.ref.value, width, height));
        }
    }

    reset() {
        if (this._supported) this._init();
    }

    private _init() {
        const { extensions: { drawBuffers } } = this.webgl;
        for (let i = 0; i < 2; i++) {
            // depth
            this.depthFramebuffers[i].bind();
            drawBuffers!.drawBuffers([
                drawBuffers!.COLOR_ATTACHMENT0,
                drawBuffers!.COLOR_ATTACHMENT1,
                drawBuffers!.COLOR_ATTACHMENT2
            ]);

            this.colorFrontTextures[i].attachFramebuffer(this.depthFramebuffers[i], 'color0');
            this.colorBackTextures[i].attachFramebuffer(this.depthFramebuffers[i], 'color1');
            this.depthTextures[i].attachFramebuffer(this.depthFramebuffers[i], 'color2');

            // color
            this.colorFramebuffers[i].bind();
            drawBuffers!.drawBuffers([
                drawBuffers!.COLOR_ATTACHMENT0,
                drawBuffers!.COLOR_ATTACHMENT1
            ]);

            this.colorFrontTextures[i].attachFramebuffer(this.colorFramebuffers[i], 'color0');
            this.colorBackTextures[i].attachFramebuffer(this.colorFramebuffers[i], 'color1');
        }
    }

    static isSupported(webgl: WebGLContext) {
        const { extensions: { drawBuffers, textureFloat, colorBufferFloat, depthTexture, blendMinMax } } = webgl;
        if (!textureFloat || !colorBufferFloat || !depthTexture || !drawBuffers || !blendMinMax) {
            if (isDebugMode) {
                const missing: string[] = [];
                if (!textureFloat) missing.push('textureFloat');
                if (!colorBufferFloat) missing.push('colorBufferFloat');
                if (!depthTexture) missing.push('depthTexture');
                if (!drawBuffers) missing.push('drawBuffers');
                if (!blendMinMax) missing.push('blendMinMax');
                console.log(`Missing "${missing.join('", "')}" extensions required for "dpoit"`);
            }
            return false;
        } else {
            return true;
        }
    }

    constructor(private webgl: WebGLContext, width: number, height: number) {
        if (!DpoitPass.isSupported(webgl)) return;

        const { resources, extensions: { colorBufferHalfFloat, textureHalfFloat } } = webgl;

        // textures

        if (isWebGL2(webgl.gl)) {
            this.depthTextures = [
                resources.texture('image-float32', 'rg', 'float', 'nearest'),
                resources.texture('image-float32', 'rg', 'float', 'nearest')
            ];

            this.colorFrontTextures = colorBufferHalfFloat && textureHalfFloat ? [
                resources.texture('image-float16', 'rgba', 'fp16', 'nearest'),
                resources.texture('image-float16', 'rgba', 'fp16', 'nearest')
            ] : [
                resources.texture('image-float32', 'rgba', 'float', 'nearest'),
                resources.texture('image-float32', 'rgba', 'float', 'nearest')
            ];

            this.colorBackTextures = colorBufferHalfFloat && textureHalfFloat ? [
                resources.texture('image-float16', 'rgba', 'fp16', 'nearest'),
                resources.texture('image-float16', 'rgba', 'fp16', 'nearest')
            ] : [
                resources.texture('image-float32', 'rgba', 'float', 'nearest'),
                resources.texture('image-float32', 'rgba', 'float', 'nearest')
            ];
        } else {
            // webgl1 requires consistent bit plane counts

            this.depthTextures = [
                resources.texture('image-float32', 'rgba', 'float', 'nearest'),
                resources.texture('image-float32', 'rgba', 'float', 'nearest')
            ];

            this.colorFrontTextures = [
                resources.texture('image-float32', 'rgba', 'float', 'nearest'),
                resources.texture('image-float32', 'rgba', 'float', 'nearest')
            ];

            this.colorBackTextures = [
                resources.texture('image-float32', 'rgba', 'float', 'nearest'),
                resources.texture('image-float32', 'rgba', 'float', 'nearest')
            ];
        }

        this.depthTextures[0].define(width, height);
        this.depthTextures[1].define(width, height);

        this.colorFrontTextures[0].define(width, height);
        this.colorFrontTextures[1].define(width, height);

        this.colorBackTextures[0].define(width, height);
        this.colorBackTextures[1].define(width, height);

        // framebuffers

        this.depthFramebuffers = [resources.framebuffer(), resources.framebuffer()];
        this.colorFramebuffers = [resources.framebuffer(), resources.framebuffer()];

        // renderables

        this.blendBackRenderable = getBlendBackDpoitRenderable(webgl, this.colorBackTextures[0]);
        this.renderable = getEvaluateDpoitRenderable(webgl, this.colorFrontTextures[0]);

        this._supported = true;
        this._init();
    }
}
