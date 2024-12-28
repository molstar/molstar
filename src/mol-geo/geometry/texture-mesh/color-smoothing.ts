/**
 * Copyright (c) 2021-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util';
import { createComputeRenderable, ComputeRenderable } from '../../../mol-gl/renderable';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { isNullTexture, Texture } from '../../../mol-gl/webgl/texture';
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
import { isWebGL2 } from '../../../mol-gl/webgl/compat';
import { TextureMeshValues } from '../../../mol-gl/renderable/texture-mesh';
import { isTimingMode } from '../../../mol-util/debug';

export const ColorAccumulateSchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),
    stride: ValueSpec('number'),

    uGroupCount: UniformSpec('i', 'material'),

    aTransform: AttributeSpec('float32', 16, 1),
    aInstance: AttributeSpec('float32', 1, 1),
    aSample: AttributeSpec('float32', 1, 0),

    uGeoTexDim: UniformSpec('v2', 'material'),
    tPosition: TextureSpec('texture', 'rgba', 'float', 'nearest', 'material'),
    tGroup: TextureSpec('texture', 'rgba', 'float', 'nearest', 'material'),

    uColorTexDim: UniformSpec('v2', 'material'),
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest', 'material'),
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
const ColorCountName = 'color-count';

interface AccumulateInput {
    vertexCount: number
    instanceCount: number
    groupCount: number
    transformBuffer: Float32Array
    instanceBuffer: Float32Array
    positionTexture: Texture
    groupTexture: Texture
    colorData: Texture
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

        const sampleCount = Math.round(input.vertexCount / stride);
        if (sampleCount > v.drawCount.ref.value || stride !== v.stride.ref.value) {
            ValueCell.update(v.aSample, getSampleBuffer(sampleCount, stride));
        }

        ValueCell.updateIfChanged(v.drawCount, sampleCount);
        ValueCell.updateIfChanged(v.instanceCount, input.instanceCount);
        ValueCell.updateIfChanged(v.stride, stride);

        ValueCell.updateIfChanged(v.uGroupCount, input.groupCount);

        ValueCell.update(v.aTransform, input.transformBuffer);
        ValueCell.update(v.aInstance, input.instanceBuffer);

        ValueCell.update(v.uGeoTexDim, Vec2.set(v.uGeoTexDim.ref.value, input.positionTexture.getWidth(), input.positionTexture.getHeight()));
        ValueCell.update(v.tPosition, input.positionTexture);
        ValueCell.update(v.tGroup, input.groupTexture);

        ValueCell.update(v.uColorTexDim, Vec2.set(v.uColorTexDim.ref.value, input.colorData.getWidth(), input.colorData.getHeight()));
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
    const sampleCount = Math.round(input.vertexCount / stride);

    const values: ColorAccumulateValues = {
        drawCount: ValueCell.create(sampleCount),
        instanceCount: ValueCell.create(input.instanceCount),
        stride: ValueCell.create(stride),

        uGroupCount: ValueCell.create(input.groupCount),

        aTransform: ValueCell.create(input.transformBuffer),
        aInstance: ValueCell.create(input.instanceBuffer),
        aSample: ValueCell.create(getSampleBuffer(sampleCount, stride)),

        uGeoTexDim: ValueCell.create(Vec2.create(input.positionTexture.getWidth(), input.positionTexture.getHeight())),
        tPosition: ValueCell.create(input.positionTexture),
        tGroup: ValueCell.create(input.groupTexture),

        uColorTexDim: ValueCell.create(Vec2.create(input.colorData.getWidth(), input.colorData.getHeight())),
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
    const shaderCode = ShaderCode('accumulate', accumulate_vert, accumulate_frag, { drawBuffers: 'required' });
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
    tCount: TextureSpec('texture', 'alpha', 'float', 'nearest'),
    uTexSize: UniformSpec('v2'),

};
type ColorNormalizeValues = Values<typeof ColorNormalizeSchema>
const ColorNormalizeName = 'color-normalize';

function getNormalizeRenderable(ctx: WebGLContext, color: Texture, count: Texture): ComputeRenderable<ColorNormalizeValues> {
    if (ctx.namedComputeRenderables[ColorNormalizeName]) {
        const v = ctx.namedComputeRenderables[ColorNormalizeName].values as ColorNormalizeValues;

        ValueCell.update(v.tColor, color);
        ValueCell.update(v.tCount, count);
        ValueCell.update(v.uTexSize, Vec2.set(v.uTexSize.ref.value, color.getWidth(), color.getHeight()));

        ctx.namedComputeRenderables[ColorNormalizeName].update();
    } else {
        ctx.namedComputeRenderables[ColorNormalizeName] = createColorNormalizeRenderable(ctx, color, count);
    }
    return ctx.namedComputeRenderables[ColorNormalizeName];
}

function createColorNormalizeRenderable(ctx: WebGLContext, color: Texture, count: Texture) {
    const values: ColorNormalizeValues = {
        ...QuadValues,
        tColor: ValueCell.create(color),
        tCount: ValueCell.create(count),
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

export function calcTextureMeshColorSmoothing(input: ColorSmoothingInput, resolution: number, stride: number, webgl: WebGLContext, texture?: Texture) {
    const { drawBuffers } = webgl.extensions;
    if (!drawBuffers) throw new Error('need WebGL draw buffers');

    if (isTimingMode) webgl.timer.mark('calcTextureMeshColorSmoothing');
    const { gl, resources, state, extensions: { colorBufferHalfFloat, textureHalfFloat } } = webgl;

    const isInstanceType = input.colorType.endsWith('Instance');
    const box = Box3D.fromSphere3D(Box3D(), isInstanceType ? input.boundingSphere : input.invariantBoundingSphere);
    const pad = 1 + resolution;
    const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));

    const scaleFactor = 1 / resolution;
    const scaledBox = Box3D.scale(Box3D(), expandedBox, scaleFactor);
    const gridDim = Box3D.size(Vec3(), scaledBox);
    Vec3.ceil(gridDim, gridDim);
    Vec3.add(gridDim, gridDim, Vec3.create(2, 2, 2));
    const { min } = expandedBox;

    const [dx, dy, dz] = gridDim;
    const { texDimX: width, texDimY: height, texCols } = getTexture2dSize(gridDim);
    // console.log({ width, height, texCols, gridDim, resolution });

    if (!webgl.namedFramebuffers[ColorAccumulateName]) {
        webgl.namedFramebuffers[ColorAccumulateName] = webgl.resources.framebuffer();
    }
    const framebuffer = webgl.namedFramebuffers[ColorAccumulateName];

    if (isWebGL2(gl)) {
        if (!webgl.namedTextures[ColorAccumulateName]) {
            webgl.namedTextures[ColorAccumulateName] = colorBufferHalfFloat && textureHalfFloat
                ? resources.texture('image-float16', 'rgba', 'fp16', 'nearest')
                : resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }

        if (!webgl.namedTextures[ColorCountName]) {
            webgl.namedTextures[ColorCountName] = resources.texture('image-float32', 'alpha', 'float', 'nearest');
        }
    } else {
        // webgl1 requires consistent bit plane counts
        // this is quite wasteful but good enough for medium size meshes

        if (!webgl.namedTextures[ColorAccumulateName]) {
            webgl.namedTextures[ColorAccumulateName] = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }

        if (!webgl.namedTextures[ColorCountName]) {
            webgl.namedTextures[ColorCountName] = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }
    }

    const accumulateTexture = webgl.namedTextures[ColorAccumulateName];
    const countTexture = webgl.namedTextures[ColorCountName];

    accumulateTexture.define(width, height);
    countTexture.define(width, height);

    accumulateTexture.attachFramebuffer(framebuffer, 0);
    countTexture.attachFramebuffer(framebuffer, 1);

    const accumulateRenderable = getAccumulateRenderable(webgl, input, expandedBox, resolution, stride);
    state.currentRenderItemId = -1;

    framebuffer.bind();
    drawBuffers.drawBuffers([
        drawBuffers.COLOR_ATTACHMENT0,
        drawBuffers.COLOR_ATTACHMENT1,
    ]);

    const { uCurrentSlice, uCurrentX, uCurrentY } = accumulateRenderable.values;

    if (isTimingMode) webgl.timer.mark('ColorAccumulate.render');
    setAccumulateDefaults(webgl);
    state.viewport(0, 0, width, height);
    state.scissor(0, 0, width, height);
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
        state.viewport(currX, currY, dx, dy);
        state.scissor(currX, currY, dx, dy);
        accumulateRenderable.render();
        ++currCol;
        currX += dx;
    }

    accumulateTexture.detachFramebuffer(framebuffer, 0);
    countTexture.detachFramebuffer(framebuffer, 1);
    drawBuffers.drawBuffers([gl.COLOR_ATTACHMENT0, gl.NONE]);
    if (isTimingMode) webgl.timer.markEnd('ColorAccumulate.render');

    // const accImage = new Float32Array(width * height * 4);
    // accumulateTexture.attachFramebuffer(framebuffer, 0);
    // webgl.readPixels(0, 0, width, height, accImage);
    // console.log(accImage);
    // printTextureImage({ array: accImage, width, height }, { scale: 1 });

    // const cntImage = new Float32Array(width * height * 4);
    // countTexture.attachFramebuffer(framebuffer, 0);
    // webgl.readPixels(0, 0, width, height, cntImage);
    // console.log(cntImage);
    // printTextureImage({ array: cntImage, width, height }, { scale: 1 });

    // normalize

    if (isTimingMode) webgl.timer.mark('ColorNormalize.render');
    if (!texture || isNullTexture(texture)) {
        texture = resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    }
    texture.define(width, height);

    const normalizeRenderable = getNormalizeRenderable(webgl, accumulateTexture, countTexture);
    state.currentRenderItemId = -1;

    setNormalizeDefaults(webgl);
    texture.attachFramebuffer(framebuffer, 0);
    state.viewport(0, 0, width, height);
    state.scissor(0, 0, width, height);
    gl.clear(gl.COLOR_BUFFER_BIT);
    normalizeRenderable.render();
    if (isTimingMode) webgl.timer.markEnd('ColorNormalize.render');

    // const normImage = new Uint8Array(width * height * 4);
    // texture.attachFramebuffer(framebuffer, 0);
    // webgl.readPixels(0, 0, width, height, normImage);
    // console.log(normImage);
    // printTextureImage({ array: normImage, width, height }, { scale: 1 });

    const gridTransform = Vec4.create(min[0], min[1], min[2], scaleFactor);
    const type = isInstanceType ? 'volumeInstance' : 'volume';
    if (isTimingMode) webgl.timer.markEnd('calcTextureMeshColorSmoothing');

    // printTextureImage(readTexture(webgl, texture), { scale: 0.75 });

    return { texture, gridDim, gridTexDim: Vec2.create(width, height), gridTransform, type };
}

//

const ColorSmoothingRgbName = 'color-smoothing-rgb';
const ColorSmoothingRgbaName = 'color-smoothing-rgba';
const ColorSmoothingAlphaName = 'color-smoothing-alpha';

function isSupportedColorType(x: string): x is 'group' | 'groupInstance' {
    return x === 'group' || x === 'groupInstance';
}

export function applyTextureMeshColorSmoothing(values: TextureMeshValues, resolution: number, stride: number, webgl: WebGLContext, colorTexture?: Texture) {
    if (!isSupportedColorType(values.dColorType.ref.value)) return;

    stride *= 3; // triple because TextureMesh is never indexed (no elements buffer)

    if (!webgl.namedTextures[ColorSmoothingRgbName]) {
        webgl.namedTextures[ColorSmoothingRgbName] = webgl.resources.texture('image-uint8', 'rgb', 'ubyte', 'nearest');
    }
    const colorData = webgl.namedTextures[ColorSmoothingRgbName];
    colorData.load(values.tColor.ref.value);

    const smoothingData = calcTextureMeshColorSmoothing({
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        colorData,
        colorType: values.dColorType.ref.value,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
    }, resolution, stride, webgl, colorTexture);

    ValueCell.updateIfChanged(values.dColorType, smoothingData.type);
    ValueCell.update(values.tColorGrid, smoothingData.texture);
    ValueCell.update(values.uColorTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
    ValueCell.update(values.uColorGridTransform, smoothingData.gridTransform);
}

function isSupportedOverpaintType(x: string): x is 'groupInstance' {
    return x === 'groupInstance';
}

export function applyTextureMeshOverpaintSmoothing(values: TextureMeshValues, resolution: number, stride: number, webgl: WebGLContext, colorTexture?: Texture) {
    if (!isSupportedOverpaintType(values.dOverpaintType.ref.value)) return;

    stride *= 3; // triple because TextureMesh is never indexed (no elements buffer)

    if (!webgl.namedTextures[ColorSmoothingRgbaName]) {
        webgl.namedTextures[ColorSmoothingRgbaName] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }
    const colorData = webgl.namedTextures[ColorSmoothingRgbaName];
    colorData.load(values.tOverpaint.ref.value);

    const smoothingData = calcTextureMeshColorSmoothing({
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        colorData,
        colorType: values.dOverpaintType.ref.value,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
    }, resolution, stride, webgl, colorTexture);

    ValueCell.updateIfChanged(values.dOverpaintType, smoothingData.type);
    ValueCell.update(values.tOverpaintGrid, smoothingData.texture);
    ValueCell.update(values.uOverpaintTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uOverpaintGridDim, smoothingData.gridDim);
    ValueCell.update(values.uOverpaintGridTransform, smoothingData.gridTransform);
}

function isSupportedTransparencyType(x: string): x is 'groupInstance' {
    return x === 'groupInstance';
}

export function applyTextureMeshTransparencySmoothing(values: TextureMeshValues, resolution: number, stride: number, webgl: WebGLContext, colorTexture?: Texture) {
    if (!isSupportedTransparencyType(values.dTransparencyType.ref.value)) return;

    stride *= 3; // triple because TextureMesh is never indexed (no elements buffer)

    if (!webgl.namedTextures[ColorSmoothingAlphaName]) {
        webgl.namedTextures[ColorSmoothingAlphaName] = webgl.resources.texture('image-uint8', 'alpha', 'ubyte', 'nearest');
    }
    const colorData = webgl.namedTextures[ColorSmoothingAlphaName];
    colorData.load(values.tTransparency.ref.value);

    const smoothingData = calcTextureMeshColorSmoothing({
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        colorData,
        colorType: values.dTransparencyType.ref.value,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
    }, resolution, stride, webgl, colorTexture);

    ValueCell.updateIfChanged(values.dTransparencyType, smoothingData.type);
    ValueCell.update(values.tTransparencyGrid, smoothingData.texture);
    ValueCell.update(values.uTransparencyTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uTransparencyGridDim, smoothingData.gridDim);
    ValueCell.update(values.uTransparencyGridTransform, smoothingData.gridTransform);
}

function isSupportedEmissiveType(x: string): x is 'groupInstance' {
    return x === 'groupInstance';
}

export function applyTextureMeshEmissiveSmoothing(values: TextureMeshValues, resolution: number, stride: number, webgl: WebGLContext, colorTexture?: Texture) {
    if (!isSupportedEmissiveType(values.dEmissiveType.ref.value)) return;

    stride *= 3; // triple because TextureMesh is never indexed (no elements buffer)

    if (!webgl.namedTextures[ColorSmoothingAlphaName]) {
        webgl.namedTextures[ColorSmoothingAlphaName] = webgl.resources.texture('image-uint8', 'alpha', 'ubyte', 'nearest');
    }
    const colorData = webgl.namedTextures[ColorSmoothingAlphaName];
    colorData.load(values.tEmissive.ref.value);

    const smoothingData = calcTextureMeshColorSmoothing({
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        colorData,
        colorType: values.dEmissiveType.ref.value,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
    }, resolution, stride, webgl, colorTexture);

    ValueCell.updateIfChanged(values.dEmissiveType, smoothingData.type);
    ValueCell.update(values.tEmissiveGrid, smoothingData.texture);
    ValueCell.update(values.uEmissiveTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uEmissiveGridDim, smoothingData.gridDim);
    ValueCell.update(values.uEmissiveGridTransform, smoothingData.gridTransform);
}

function isSupportedSubstanceType(x: string): x is 'groupInstance' {
    return x === 'groupInstance';
}

export function applyTextureMeshSubstanceSmoothing(values: TextureMeshValues, resolution: number, stride: number, webgl: WebGLContext, colorTexture?: Texture) {
    if (!isSupportedSubstanceType(values.dSubstanceType.ref.value)) return;

    stride *= 3; // triple because TextureMesh is never indexed (no elements buffer)

    if (!webgl.namedTextures[ColorSmoothingRgbaName]) {
        webgl.namedTextures[ColorSmoothingRgbaName] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }
    const colorData = webgl.namedTextures[ColorSmoothingRgbaName];
    colorData.load(values.tSubstance.ref.value);

    const smoothingData = calcTextureMeshColorSmoothing({
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        colorData,
        colorType: values.dSubstanceType.ref.value,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
    }, resolution, stride, webgl, colorTexture);

    ValueCell.updateIfChanged(values.dSubstanceType, smoothingData.type);
    ValueCell.update(values.tSubstanceGrid, smoothingData.texture);
    ValueCell.update(values.uSubstanceTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uSubstanceGridDim, smoothingData.gridDim);
    ValueCell.update(values.uSubstanceGridTransform, smoothingData.gridTransform);
}
