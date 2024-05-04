/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Partially adapted from three.js, The MIT License, Copyright © 2010-2024 three.js authors
 */

import { CopyRenderable, QuadSchema, QuadValues, createCopyRenderable } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture, createNullTexture } from '../../mol-gl/webgl/texture';
import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { isTimingMode } from '../../mol-util/debug';
import { composite_frag } from '../../mol-gl/shader/bloom/composite.frag';
import { luminosity_frag } from '../../mol-gl/shader/bloom/luminosity.frag';
import { blur_frag } from '../../mol-gl/shader/bloom/blur.frag';
import { memoize1 } from '../../mol-util/memoize';
import { PostprocessingProps } from './postprocessing';

const MipCount = 5;

export const BloomParams = {
    strength: PD.Numeric(1, { min: 0, max: 3, step: 0.1 }),
    radius: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
    threshold: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { description: 'Luminosity threshold', hideIf: p => p.mode === 'emissive' }),
    mode: PD.Select('emissive', [['luminosity', 'Luminosity'], ['emissive', 'Emissive']] as const),
};
export type BloomProps = PD.Values<typeof BloomParams>

export class BloomPass {
    static isEnabled(props: PostprocessingProps) {
        return props.bloom.name === 'on';
    }

    readonly emissiveTarget: RenderTarget;

    private readonly luminosityTarget: RenderTarget;
    private readonly horizontalBlurTargets: RenderTarget[] = [];
    private readonly verticalBlurTargets: RenderTarget[] = [];
    private readonly compositeTarget: RenderTarget;

    private readonly luminosityRenderable: LuminosityRenderable;
    private readonly blurRenderable: BlurRenderable;
    private readonly compositeRenderable: CompositeRenderable;
    private readonly copyRenderable: CopyRenderable;

    constructor(private webgl: WebGLContext, width: number, height: number) {
        this.emissiveTarget = webgl.createRenderTarget(width, height, true, 'uint8', 'linear', 'rgba');

        this.luminosityTarget = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');
        this.compositeTarget = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');

        let blurWidth = Math.round(width / 2);
        let blurHeight = Math.round(height / 2);
        for (let i = 0; i < MipCount; ++i) {
            this.horizontalBlurTargets[i] = webgl.createRenderTarget(blurWidth, blurHeight, false, 'uint8', 'linear');
            this.verticalBlurTargets[i] = webgl.createRenderTarget(blurWidth, blurHeight, false, 'uint8', 'linear');
            blurWidth = Math.round(blurWidth / 2);
            blurHeight = Math.round(blurHeight / 2);
        }

        const nullTexture = createNullTexture();
        this.luminosityRenderable = getLuminosityRenderable(webgl, nullTexture, nullTexture, nullTexture);
        this.blurRenderable = getBlurRenderable(webgl, nullTexture);
        this.compositeRenderable = getCompositeRenderable(webgl, width, height, this.verticalBlurTargets[0].texture, this.verticalBlurTargets[1].texture, this.verticalBlurTargets[2].texture, this.verticalBlurTargets[3].texture, this.verticalBlurTargets[4].texture);
        this.copyRenderable = createCopyRenderable(webgl, this.compositeTarget.texture);
    }

    setSize(width: number, height: number) {
        const w = this.luminosityTarget.getWidth();
        const h = this.luminosityTarget.getHeight();

        if (width !== w || height !== h) {
            this.emissiveTarget.setSize(width, height);
            this.luminosityTarget.setSize(width, height);
            this.compositeTarget.setSize(width, height);

            let blurWidth = Math.round(width / 2);
            let blurHeight = Math.round(height / 2);
            for (let i = 0; i < MipCount; ++i) {
                this.horizontalBlurTargets[i].setSize(blurWidth, blurHeight);
                this.verticalBlurTargets[i].setSize(blurWidth, blurHeight);
                blurWidth = Math.round(blurWidth / 2);
                blurHeight = Math.round(blurHeight / 2);
            }

            ValueCell.update(this.luminosityRenderable.values.uTexSizeInv, Vec2.set(this.compositeRenderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.compositeRenderable.values.uTexSizeInv, Vec2.set(this.compositeRenderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));

            ValueCell.update(this.copyRenderable.values.uTexSize, Vec2.set(this.copyRenderable.values.uTexSize.ref.value, width, height));
        }
    }

    update(input: Texture, emissive: Texture, depth: Texture, props: BloomProps) {
        let luminosityNeedsUpdate = false;

        if (this.luminosityRenderable.values.tColor.ref.value !== input) {
            ValueCell.update(this.luminosityRenderable.values.tColor, input);
            luminosityNeedsUpdate = true;
        }

        if (this.luminosityRenderable.values.tEmissive.ref.value !== emissive) {
            ValueCell.update(this.luminosityRenderable.values.tEmissive, emissive);
            luminosityNeedsUpdate = true;
        }

        if (this.luminosityRenderable.values.tDepth.ref.value !== depth) {
            ValueCell.update(this.luminosityRenderable.values.tDepth, depth);
            luminosityNeedsUpdate = true;
        }

        if (this.luminosityRenderable.values.dMode.ref.value !== props.mode) {
            ValueCell.update(this.luminosityRenderable.values.dMode, props.mode);
            luminosityNeedsUpdate = true;
        }

        ValueCell.updateIfChanged(this.luminosityRenderable.values.uLuminosityThreshold, props.threshold);

        if (luminosityNeedsUpdate) {
            this.luminosityRenderable.update();
        }

        //

        let blurNeedsUpdate = false;

        if (this.blurRenderable.values.tInput.ref.value !== input) {
            ValueCell.update(this.blurRenderable.values.tInput, input);
            blurNeedsUpdate = true;
        }

        if (blurNeedsUpdate) {
            this.blurRenderable.update();
        }

        //

        ValueCell.update(this.compositeRenderable.values.uBloomRadius, props.radius);
        ValueCell.update(this.compositeRenderable.values.uBloomStrength, props.strength);
    }

    render(viewport: Viewport, target: RenderTarget | undefined) {
        if (isTimingMode) this.webgl.timer.mark('BloomPass.render');
        const { gl, state } = this.webgl;
        const { x, y, width, height } = viewport;
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        // printTextureImage(readTexture(this.webgl, this.luminosityRenderable.values.tEmissive.ref.value, new Uint8Array(width * height * 4)), { scale: 0.25, id: 'emissive' });

        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        // Extract Bright Areas
        this.luminosityTarget.bind();
        this.luminosityRenderable.render();

        // printTextureImage(readTexture(this.webgl, this.luminosityTarget.texture, new Uint8Array(width * height * 4)), { scale: 0.25, id: 'luminosity' });

        // Blur All the mips progressively
        for (let i = 0; i < MipCount; ++i) {
            const blurWidth = this.horizontalBlurTargets[i].getWidth();
            const blurHeight = this.horizontalBlurTargets[i].getHeight();
            state.viewport(0, 0, blurWidth, blurHeight);
            state.scissor(0, 0, blurWidth, blurHeight);

            ValueCell.update(this.blurRenderable.values.dKernelRadius, BlurKernelSizes[i]);
            ValueCell.update(this.blurRenderable.values.uGaussianCoefficients, getBlurCoefficients(BlurKernelSizes[i]));
            ValueCell.update(this.blurRenderable.values.uTexSizeInv, Vec2.set(this.blurRenderable.values.uTexSizeInv.ref.value, 1 / blurWidth, 1 / blurHeight));

            this.horizontalBlurTargets[i].bind();
            ValueCell.update(this.blurRenderable.values.tInput, i === 0 ? this.luminosityTarget.texture : this.verticalBlurTargets[i - 1].texture);
            ValueCell.update(this.blurRenderable.values.uDirection, BlurDirectionX);
            this.blurRenderable.update();
            this.blurRenderable.render();

            this.verticalBlurTargets[i].bind();
            ValueCell.update(this.blurRenderable.values.tInput, this.horizontalBlurTargets[i].texture);
            ValueCell.update(this.blurRenderable.values.uDirection, BlurDirectionY);
            this.blurRenderable.update();
            this.blurRenderable.render();

            // printTextureImage(readTexture(this.webgl, this.verticalBlurTargets[i].texture, new Uint8Array(blurWidth * blurHeight * 4)), { scale: 0.25, id: `blur-${i}` });
        }

        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        // Composite All the mips
        this.compositeTarget.bind();
        this.compositeRenderable.update();
        this.compositeRenderable.render();

        // printTextureImage(readTexture(this.webgl, this.compositeTarget.texture, new Uint8Array(width * height * 4)), { scale: 0.25, id: 'composite' });

        if (target) {
            target.bind();
        } else {
            this.webgl.unbindFramebuffer();
        }
        state.enable(gl.BLEND);
        state.blendFunc(gl.ONE, gl.ONE);
        this.copyRenderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('BloomPass.render');
    }
}

//

const LuminositySchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tEmissive: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSizeInv: UniformSpec('v2'),

    uDefaultColor: UniformSpec('v3'),
    uDefaultOpacity: UniformSpec('f'),
    uLuminosityThreshold: UniformSpec('f'),
    uSmoothWidth: UniformSpec('f'),
    dMode: DefineSpec('string', ['luminosity', 'emissive']),
};
const LuminosityShaderCode = ShaderCode('Bloom Luminosity', quad_vert, luminosity_frag);
type LuminosityRenderable = ComputeRenderable<Values<typeof LuminositySchema>>

function getLuminosityRenderable(ctx: WebGLContext, colorTexture: Texture, emissiveTexture: Texture, depthTexture: Texture): LuminosityRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof LuminositySchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        tEmissive: ValueCell.create(emissiveTexture),
        tDepth: ValueCell.create(depthTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

        uDefaultColor: ValueCell.create(Vec3()),
        uDefaultOpacity: ValueCell.create(0),
        uLuminosityThreshold: ValueCell.create(0),
        uSmoothWidth: ValueCell.create(1),
        dMode: ValueCell.create('emissive'),
    };

    const schema = { ...LuminositySchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', LuminosityShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

function _getBlurCoefficients(kernelRadius: number) {
    const coefficients: number[] = [];
    for (let i = 0; i < kernelRadius; ++i) {
        coefficients.push(0.39894 * Math.exp(-0.5 * i * i / (kernelRadius * kernelRadius)) / kernelRadius);
    }
    return coefficients;
}
const getBlurCoefficients = memoize1(_getBlurCoefficients);

const BlurKernelSizes = [3, 5, 7, 9, 11];

const BlurDirectionX = Vec2.create(1, 0);
const BlurDirectionY = Vec2.create(0, 1);

const BlurSchema = {
    ...QuadSchema,
    tInput: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),

    uDirection: UniformSpec('v2'),
    uGaussianCoefficients: UniformSpec('f[]'),
    dKernelRadius: DefineSpec('number'),
};
const BlurShaderCode = ShaderCode('Bloom Blur', quad_vert, blur_frag);
type BlurRenderable = ComputeRenderable<Values<typeof BlurSchema>>

function getBlurRenderable(ctx: WebGLContext, inputTexture: Texture): BlurRenderable {
    const width = inputTexture.getWidth();
    const height = inputTexture.getHeight();

    const values: Values<typeof BlurSchema> = {
        ...QuadValues,
        tInput: ValueCell.create(inputTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

        uDirection: ValueCell.create(Vec2()),
        uGaussianCoefficients: ValueCell.create([]),
        dKernelRadius: ValueCell.create(BlurKernelSizes[0]),
    };

    const schema = { ...BlurSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', BlurShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

const CompositeSchema = {
    ...QuadSchema,
    tBlur1: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tBlur2: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tBlur3: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tBlur4: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tBlur5: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),

    uBloomStrength: UniformSpec('f'),
    uBloomRadius: UniformSpec('f'),
    uBloomFactors: UniformSpec('f[]'),
    uBloomTints: UniformSpec('v3[]'),
};
const CompositeShaderCode = ShaderCode('Bloom Composite', quad_vert, composite_frag);
type CompositeRenderable = ComputeRenderable<Values<typeof CompositeSchema>>

function getCompositeRenderable(ctx: WebGLContext, width: number, height: number, blurTexture1: Texture, blurTexture2: Texture, blurTexture3: Texture, blurTexture4: Texture, blurTexture5: Texture): CompositeRenderable {
    const values: Values<typeof CompositeSchema> = {
        ...QuadValues,
        uTexSizeInv: ValueCell.create(Vec2.create(width, height)),

        tBlur1: ValueCell.create(blurTexture1),
        tBlur2: ValueCell.create(blurTexture2),
        tBlur3: ValueCell.create(blurTexture3),
        tBlur4: ValueCell.create(blurTexture4),
        tBlur5: ValueCell.create(blurTexture5),

        uBloomStrength: ValueCell.create(1),
        uBloomRadius: ValueCell.create(0),
        uBloomFactors: ValueCell.create([1.0, 0.8, 0.6, 0.4, 0.2]),
        uBloomTints: ValueCell.create(new Array(5 * 3).fill(1)),
    };

    const schema = { ...CompositeSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', CompositeShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}