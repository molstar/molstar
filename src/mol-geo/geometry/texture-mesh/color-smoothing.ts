/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util';
import { createComputeRenderable, ComputeRenderable } from '../../../mol-gl/renderable';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Texture } from '../../../mol-gl/webgl/texture';
import { ShaderCode } from '../../../mol-gl/shader-code';
import { createComputeRenderItem } from '../../../mol-gl/webgl/render-item';
import { ValueSpec, AttributeSpec, UniformSpec, TextureSpec, Values, DefineSpec } from '../../../mol-gl/renderable/schema';
import { quad_vert } from '../../../mol-gl/shader/quad.vert';
import { normalize_frag } from '../../../mol-gl/shader/compute/color-smoothing/normalize.frag';
import { QuadSchema, QuadValues } from '../../../mol-gl/compute/util';
import { Vec2, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { Box3D, Sphere3D } from '../../../mol-math/geometry';
import { accumulate_frag } from '../../../mol-gl/shader/compute/color-smoothing/accumulate.frag';
import { accumulate_vert } from '../../../mol-gl/shader/compute/color-smoothing/accumulate.vert';
import { TextureImage } from '../../../mol-gl/renderable/util';

export const ColorAccumulateSchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    uTotalCount: UniformSpec('i'),
    uInstanceCount: UniformSpec('i'),
    uGroupCount: UniformSpec('i'),

    aTransform: AttributeSpec('float32', 16, 1),
    aInstance: AttributeSpec('float32', 1, 1),
    aSample: AttributeSpec('float32', 1, 0),

    uGeoTexDim: UniformSpec('v2', 'buffered'),
    tPosition: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tGroup: TextureSpec('texture', 'rgba', 'float', 'nearest'),

    uColorTexDim: UniformSpec('v2'),
    tColor: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    dColorType: DefineSpec('string', ['group', 'groupInstance', 'vertex', 'vertexInstance']),

    uCurrentSlice: UniformSpec('f'),
    uCurrentX: UniformSpec('f'),
    uCurrentY: UniformSpec('f'),
    uBboxMin: UniformSpec('v3', 'material'),
    uBboxSize: UniformSpec('v3', 'material'),
    uResolution: UniformSpec('f', 'material'),
};
type ColorAccumulateValues = Values<typeof ColorAccumulateSchema>
const ColorAccumulateName = 'color-accumulate';

interface AccumulateInput {
    vertexCount: number
    instanceCount: number
    groupCount: number
    transformBuffer: Float32Array
    instanceBuffer: Float32Array
    positionTexture: Texture
    groupTexture: Texture
    colorData: TextureImage<Uint8Array>
    colorType: 'group' | 'groupInstance'
}

function getSampleBuffer(sampleCount: number, stride: number) {
    const sampleBuffer = new Float32Array(sampleCount);
    for (let i = 0; i < sampleCount; ++i) {
        sampleBuffer[i] = i * stride;
    }
    return sampleBuffer;
}

function getAccumulateRenderable(ctx: WebGLContext, input: AccumulateInput, box: Box3D, resolution: number, stride: number): ComputeRenderable<ColorAccumulateValues> {
    if (ctx.namedComputeRenderables[ColorAccumulateName]) {
        const extent = Vec3.sub(Vec3(), box.max, box.min);
        const v = ctx.namedComputeRenderables[ColorAccumulateName].values as ColorAccumulateValues;

        const sampleCount = input.vertexCount / stride;
        if (sampleCount > v.drawCount.ref.value) {
            ValueCell.update(v.aSample, getSampleBuffer(sampleCount, stride));
        }

        ValueCell.updateIfChanged(v.drawCount, sampleCount);
        ValueCell.updateIfChanged(v.instanceCount, input.instanceCount);

        ValueCell.updateIfChanged(v.uTotalCount, input.vertexCount);
        ValueCell.updateIfChanged(v.uInstanceCount, input.instanceCount);
        ValueCell.updateIfChanged(v.uGroupCount, input.groupCount);

        ValueCell.update(v.aTransform, input.transformBuffer);
        ValueCell.update(v.aInstance, input.instanceBuffer);

        ValueCell.update(v.uGeoTexDim, Vec2.set(v.uGeoTexDim.ref.value, input.positionTexture.getWidth(), input.positionTexture.getHeight()));
        ValueCell.update(v.tPosition, input.positionTexture);
        ValueCell.update(v.tGroup, input.groupTexture);

        ValueCell.update(v.uColorTexDim, Vec2.set(v.uColorTexDim.ref.value, input.colorData.width, input.colorData.height));
        ValueCell.update(v.tColor, input.colorData);
        ValueCell.updateIfChanged(v.dColorType, input.colorType);

        ValueCell.updateIfChanged(v.uCurrentSlice, 0);
        ValueCell.updateIfChanged(v.uCurrentX, 0);
        ValueCell.updateIfChanged(v.uCurrentY, 0);
        ValueCell.update(v.uBboxMin, box.min);
        ValueCell.update(v.uBboxSize, extent);
        ValueCell.updateIfChanged(v.uResolution, resolution);

        ctx.namedComputeRenderables[ColorAccumulateName].update();
    } else {
        ctx.namedComputeRenderables[ColorAccumulateName] = createAccumulateRenderable(ctx, input, box, resolution, stride);
    }
    return ctx.namedComputeRenderables[ColorAccumulateName];
}

function createAccumulateRenderable(ctx: WebGLContext, input: AccumulateInput, box: Box3D, resolution: number, stride: number) {
    const extent = Vec3.sub(Vec3(), box.max, box.min);
    const sampleCount = input.vertexCount / stride;

    const values: ColorAccumulateValues = {
        drawCount: ValueCell.create(sampleCount),
        instanceCount: ValueCell.create(input.instanceCount),

        uTotalCount: ValueCell.create(input.vertexCount),
        uInstanceCount: ValueCell.create(input.instanceCount),
        uGroupCount: ValueCell.create(input.groupCount),

        aTransform: ValueCell.create(input.transformBuffer),
        aInstance: ValueCell.create(input.instanceBuffer),
        aSample: ValueCell.create(getSampleBuffer(sampleCount, stride)),

        uGeoTexDim: ValueCell.create(Vec2.create(input.positionTexture.getWidth(), input.positionTexture.getHeight())),
        tPosition: ValueCell.create(input.positionTexture),
        tGroup: ValueCell.create(input.groupTexture),

        uColorTexDim: ValueCell.create(Vec2.create(input.colorData.width, input.colorData.height)),
        tColor: ValueCell.create(input.colorData),
        dColorType: ValueCell.create(input.colorType),

        uCurrentSlice: ValueCell.create(0),
        uCurrentX: ValueCell.create(0),
        uCurrentY: ValueCell.create(0),
        uBboxMin: ValueCell.create(box.min),
        uBboxSize: ValueCell.create(extent),
        uResolution: ValueCell.create(resolution),
    };

    const schema = { ...ColorAccumulateSchema };
    const shaderCode = ShaderCode('accumulate', accumulate_vert, accumulate_frag);
    const renderItem = createComputeRenderItem(ctx, 'points', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function setAccumulateDefaults(ctx: WebGLContext) {
    const { gl, state } = ctx;
    state.disable(gl.CULL_FACE);
    state.enable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.enable(gl.SCISSOR_TEST);
    state.depthMask(false);
    state.clearColor(0, 0, 0, 0);
    state.blendFunc(gl.ONE, gl.ONE);
    state.blendEquation(gl.FUNC_ADD);
}

//

export const ColorNormalizeSchema = {
    ...QuadSchema,

    tColor: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    uTexSize: UniformSpec('v2'),

};
type ColorNormalizeValues = Values<typeof ColorNormalizeSchema>
const ColorNormalizeName = 'color-normalize';

function getNormalizeRenderable(ctx: WebGLContext, color: Texture): ComputeRenderable<ColorNormalizeValues> {
    if (ctx.namedComputeRenderables[ColorNormalizeName]) {
        const v = ctx.namedComputeRenderables[ColorNormalizeName].values as ColorNormalizeValues;

        ValueCell.update(v.tColor, color);
        ValueCell.update(v.uTexSize, Vec2.set(v.uTexSize.ref.value, color.getWidth(), color.getHeight()));

        ctx.namedComputeRenderables[ColorNormalizeName].update();
    } else {
        ctx.namedComputeRenderables[ColorNormalizeName] = createColorNormalizeRenderable(ctx, color);
    }
    return ctx.namedComputeRenderables[ColorNormalizeName];
}

function createColorNormalizeRenderable(ctx: WebGLContext, color: Texture) {
    const values: ColorNormalizeValues = {
        ...QuadValues,
        tColor: ValueCell.create(color),
        uTexSize: ValueCell.create(Vec2.create(color.getWidth(), color.getHeight())),
    };

    const schema = { ...ColorNormalizeSchema };
    const shaderCode = ShaderCode('normalize', quad_vert, normalize_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function setNormalizeDefaults(ctx: WebGLContext) {
    const { gl, state } = ctx;
    state.disable(gl.CULL_FACE);
    state.enable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.enable(gl.SCISSOR_TEST);
    state.depthMask(false);
    state.clearColor(0, 0, 0, 0);
    state.blendFunc(gl.ONE, gl.ONE);
    state.blendEquation(gl.FUNC_ADD);
}

//

function getTexture2dSize(gridDim: Vec3) {
    const area = gridDim[0] * gridDim[1] * gridDim[2];
    const squareDim = Math.sqrt(area);
    const powerOfTwoSize = Math.pow(2, Math.ceil(Math.log(squareDim) / Math.log(2)));

    let texDimX = 0;
    let texDimY = gridDim[1];
    let texRows = 1;
    let texCols = gridDim[2];
    if (powerOfTwoSize < gridDim[0] * gridDim[2]) {
        texCols = Math.floor(powerOfTwoSize / gridDim[0]);
        texRows = Math.ceil(gridDim[2] / texCols);
        texDimX = texCols * gridDim[0];
        texDimY *= texRows;
    } else {
        texDimX = gridDim[0] * gridDim[2];
    }
    // console.log(texDimX, texDimY, texDimY < powerOfTwoSize ? powerOfTwoSize : powerOfTwoSize * 2);
    return { texDimX, texDimY, texRows, texCols, powerOfTwoSize: texDimY < powerOfTwoSize ? powerOfTwoSize : powerOfTwoSize * 2 };
}

interface ColorSmoothingInput extends AccumulateInput {
    boundingSphere: Sphere3D
    invariantBoundingSphere: Sphere3D
}

export function calcTextureMeshColorSmoothing(webgl: WebGLContext, input: ColorSmoothingInput, resolution: number, stride: number, texture?: Texture) {
    const { gl, resources, state, extensions: { colorBufferHalfFloat, textureHalfFloat } } = webgl;

    const isInstanceType = input.colorType.endsWith('Instance');
    const box = Box3D.fromSphere3D(Box3D(), isInstanceType ? input.boundingSphere : input.invariantBoundingSphere);

    const scaleFactor = 1 / resolution;
    const scaledBox = Box3D.scale(Box3D(), box, scaleFactor);
    const dim = Box3D.size(Vec3(), scaledBox);
    Vec3.ceil(dim, dim);
    Vec3.add(dim, dim, Vec3.create(2, 2, 2));
    const { min } = box;

    const [ dx, dy, dz ] = dim;
    const { texDimX: width, texDimY: height, texCols } = getTexture2dSize(dim);
    // console.log({ width, height, texCols, dim, resolution });

    if (!webgl.namedTextures[ColorAccumulateName]) {
        webgl.namedTextures[ColorAccumulateName] = colorBufferHalfFloat && textureHalfFloat
            ? resources.texture('image-float16', 'rgba', 'fp16', 'nearest')
            : resources.texture('image-float32', 'rgba', 'float', 'nearest');
    }
    const accumulateTexture = webgl.namedTextures[ColorAccumulateName];
    accumulateTexture.define(width, height);

    const accumulateRenderable = getAccumulateRenderable(webgl, input, box, resolution, stride);

    //

    const { uCurrentSlice, uCurrentX, uCurrentY } = accumulateRenderable.values;

    if (!webgl.namedFramebuffers[ColorAccumulateName]) {
        webgl.namedFramebuffers[ColorAccumulateName] = webgl.resources.framebuffer();
    }
    const framebuffer = webgl.namedFramebuffers[ColorAccumulateName];
    framebuffer.bind();

    setAccumulateDefaults(webgl);
    state.currentRenderItemId = -1;
    accumulateTexture.attachFramebuffer(framebuffer, 0);
    gl.viewport(0, 0, width, height);
    gl.scissor(0, 0, width, height);
    gl.clear(gl.COLOR_BUFFER_BIT);
    ValueCell.update(uCurrentY, 0);
    let currCol = 0;
    let currY = 0;
    let currX = 0;
    for (let i = 0; i < dz; ++i) {
        if (currCol >= texCols) {
            currCol -= texCols;
            currY += dy;
            currX = 0;
            ValueCell.update(uCurrentY, currY);
        }
        // console.log({ i, currX, currY });
        ValueCell.update(uCurrentX, currX);
        ValueCell.update(uCurrentSlice, i);
        gl.viewport(currX, currY, dx, dy);
        gl.scissor(currX, currY, dx, dy);
        accumulateRenderable.render();
        ++currCol;
        currX += dx;
    }

    // const accImage = new Float32Array(width * height * 4);
    // accumulateTexture.attachFramebuffer(framebuffer, 0);
    // webgl.readPixels(0, 0, width, height, accImage);
    // console.log(accImage);
    // printTextureImage({ array: accImage, width, height }, 1 / 4);

    // normalize

    if (!texture) texture = resources.texture('image-uint8', 'rgb', 'ubyte', 'linear');
    texture.define(width, height);

    const normalizeRenderable = getNormalizeRenderable(webgl, accumulateTexture);

    setNormalizeDefaults(webgl);
    state.currentRenderItemId = -1;
    texture.attachFramebuffer(framebuffer, 0);
    gl.viewport(0, 0, width, height);
    gl.scissor(0, 0, width, height);
    gl.clear(gl.COLOR_BUFFER_BIT);
    normalizeRenderable.render();

    // const normImage = new Uint8Array(width * height * 4);
    // texture.attachFramebuffer(framebuffer, 0);
    // webgl.readPixels(0, 0, width, height, normImage);
    // console.log(normImage);
    // printTextureImage({ array: normImage, width, height }, 1 / 4);

    const transform = Vec4.create(min[0], min[1], min[2], scaleFactor);
    const type = isInstanceType ? 'volumeInstance' : 'volume';

    return { texture, gridDim: dim, gridTexDim: Vec2.create(width, height), transform, type };
}
