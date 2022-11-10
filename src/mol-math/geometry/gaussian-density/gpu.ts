/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import { PositionData } from '../common';
import { Box3D } from '../../geometry';
import { GaussianDensityProps, GaussianDensityData, GaussianDensityTextureData } from '../gaussian-density';
import { OrderedSet } from '../../../mol-data/int';
import { Vec3, Tensor, Mat4, Vec2 } from '../../linear-algebra';
import { ValueCell } from '../../../mol-util';
import { createComputeRenderable, ComputeRenderable } from '../../../mol-gl/renderable';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Texture, TextureFilter, TextureFormat, TextureKind, TextureType } from '../../../mol-gl/webgl/texture';
import { unpackRGBToInt } from '../../../mol-util/number-packing';
import { ShaderCode } from '../../../mol-gl/shader-code';
import { createComputeRenderItem } from '../../../mol-gl/webgl/render-item';
import { ValueSpec, AttributeSpec, UniformSpec, TextureSpec, DefineSpec, Values } from '../../../mol-gl/renderable/schema';
import { gaussianDensity_vert } from '../../../mol-gl/shader/gaussian-density.vert';
import { gaussianDensity_frag } from '../../../mol-gl/shader/gaussian-density.frag';
import { Framebuffer } from '../../../mol-gl/webgl/framebuffer';
import { isTimingMode } from '../../../mol-util/debug';

const GaussianDensitySchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    aRadius: AttributeSpec('float32', 1, 0),
    aPosition: AttributeSpec('float32', 3, 0),
    aGroup: AttributeSpec('float32', 1, 0),

    uCurrentSlice: UniformSpec('f'),
    uCurrentX: UniformSpec('f'),
    uCurrentY: UniformSpec('f'),
    uBboxMin: UniformSpec('v3', 'material'),
    uBboxSize: UniformSpec('v3', 'material'),
    uGridDim: UniformSpec('v3', 'material'),
    uGridTexDim: UniformSpec('v3', 'material'),
    uGridTexScale: UniformSpec('v2', 'material'),
    uAlpha: UniformSpec('f', 'material'),
    uResolution: UniformSpec('f', 'material'),
    uRadiusFactorInv: UniformSpec('f', 'material'),
    tMinDistanceTex: TextureSpec('texture', 'rgba', 'float', 'nearest', 'material'),

    dGridTexType: DefineSpec('string', ['2d', '3d']),
    dCalcType: DefineSpec('string', ['density', 'minDistance', 'groupId']),
};
type GaussianDensityValues = Values<typeof GaussianDensitySchema>
type GaussianDensityRenderable = ComputeRenderable<GaussianDensityValues>
const GaussianDensityName = 'gaussian-density';

function getFramebuffer(webgl: WebGLContext): Framebuffer {
    if (!webgl.namedFramebuffers[GaussianDensityName]) {
        webgl.namedFramebuffers[GaussianDensityName] = webgl.resources.framebuffer();
    }
    return webgl.namedFramebuffers[GaussianDensityName];
}

function getTexture(name: string, webgl: WebGLContext, kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter): Texture {
    const _name = `${GaussianDensityName}-${name}`;
    if (!webgl.namedTextures[_name]) {
        webgl.namedTextures[_name] = webgl.resources.texture(kind, format, type, filter);
    }
    return webgl.namedTextures[_name];
}

export function GaussianDensityGPU(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, webgl: WebGLContext): GaussianDensityData {
    // always use texture2d when the gaussian density needs to be downloaded from the GPU,
    // it's faster than texture3d
    // console.time('GaussianDensityTexture2d')
    const tmpTexture = getTexture('tmp', webgl, 'image-uint8', 'rgba', 'ubyte', 'linear');
    const { scale, bbox, texture, gridDim, gridTexDim, radiusFactor, resolution, maxRadius } = calcGaussianDensityTexture2d(webgl, position, box, radius, false, props, tmpTexture);
    // webgl.waitForGpuCommandsCompleteSync()
    // console.timeEnd('GaussianDensityTexture2d')
    const { field, idField } = fieldFromTexture2d(webgl, texture, gridDim, gridTexDim);

    return { field, idField, transform: getTransform(scale, bbox), radiusFactor, resolution, maxRadius };
}

export function GaussianDensityTexture(webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, oldTexture?: Texture): GaussianDensityTextureData {
    return webgl.isWebGL2 ?
        GaussianDensityTexture3d(webgl, position, box, radius, props, oldTexture) :
        GaussianDensityTexture2d(webgl, position, box, radius, false, props, oldTexture);
}

export function GaussianDensityTexture2d(webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, powerOfTwo: boolean, props: GaussianDensityProps, oldTexture?: Texture): GaussianDensityTextureData {
    if (isTimingMode) webgl.timer.mark('GaussianDensityTexture2d');
    const data = calcGaussianDensityTexture2d(webgl, position, box, radius, powerOfTwo, props, oldTexture);
    if (isTimingMode) webgl.timer.markEnd('GaussianDensityTexture2d');
    return finalizeGaussianDensityTexture(data);
}

export function GaussianDensityTexture3d(webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, oldTexture?: Texture): GaussianDensityTextureData {
    if (isTimingMode) webgl.timer.mark('GaussianDensityTexture3d');
    const data = calcGaussianDensityTexture3d(webgl, position, box, radius, props, oldTexture);
    if (isTimingMode) webgl.timer.markEnd('GaussianDensityTexture3d');
    return finalizeGaussianDensityTexture(data);
}

function finalizeGaussianDensityTexture({ texture, scale, bbox, gridDim, gridTexDim, gridTexScale, radiusFactor, resolution, maxRadius }: _GaussianDensityTextureData): GaussianDensityTextureData {
    return { transform: getTransform(scale, bbox), texture, bbox, gridDim, gridTexDim, gridTexScale, radiusFactor, resolution, maxRadius };
}

function getTransform(scale: Vec3, bbox: Box3D) {
    const transform = Mat4.identity();
    Mat4.fromScaling(transform, scale);
    Mat4.setTranslation(transform, bbox.min);
    return transform;
}

//

type _GaussianDensityTextureData = {
    texture: Texture,
    scale: Vec3,
    bbox: Box3D,
    gridDim: Vec3,
    gridTexDim: Vec3
    gridTexScale: Vec2
    radiusFactor: number
    resolution: number
    maxRadius: number
}

function calcGaussianDensityTexture2d(webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, powerOfTwo: boolean, props: GaussianDensityProps, texture?: Texture): _GaussianDensityTextureData {
    // console.log('2d');
    const { gl, resources, state, extensions: { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat, blendMinMax } } = webgl;
    const { smoothness, resolution } = props;

    const { drawCount, positions, radii, groups, scale, expandedBox, dim, maxRadius } = prepareGaussianDensityData(position, box, radius, props);
    const [dx, dy, dz] = dim;
    const { texDimX, texDimY, texCols, powerOfTwoSize } = getTexture2dSize(dim);
    // console.log({ texDimX, texDimY, texCols, powerOfTwoSize, dim });
    const gridTexDim = Vec3.create(texDimX, texDimY, 0);
    const gridTexScale = Vec2.create(texDimX / powerOfTwoSize, texDimY / powerOfTwoSize);
    const radiusFactor = maxRadius * 2;

    const width = powerOfTwo ? powerOfTwoSize : texDimX;
    const height = powerOfTwo ? powerOfTwoSize : texDimY;

    const minDistTex = getTexture('min-dist-2d', webgl, 'image-uint8', 'rgba', 'ubyte', 'nearest');
    minDistTex.define(width, height);

    const renderable = getGaussianDensityRenderable(webgl, drawCount, positions, radii, groups, minDistTex, expandedBox, dim, gridTexDim, gridTexScale, smoothness, resolution, radiusFactor);

    //

    const { uCurrentSlice, uCurrentX, uCurrentY } = renderable.values;

    const framebuffer = getFramebuffer(webgl);
    framebuffer.bind();
    setRenderingDefaults(webgl);

    if (!texture) texture = colorBufferHalfFloat && textureHalfFloat
        ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
        : colorBufferFloat && textureFloat
            ? resources.texture('image-float32', 'rgba', 'float', 'linear')
            : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    texture.define(width, height);

    // console.log(renderable)

    function render(fbTex: Texture, clear: boolean) {
        state.currentRenderItemId = -1;
        fbTex.attachFramebuffer(framebuffer, 0);
        if (clear) {
            state.viewport(0, 0, width, height);
            state.scissor(0, 0, width, height);
            gl.clear(gl.COLOR_BUFFER_BIT);
        }
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
            renderable.render();
            ++currCol;
            currX += dx;
        }
        gl.flush();
    }

    setupDensityRendering(webgl, renderable);
    render(texture, true);

    if (blendMinMax) {
        setupMinDistanceRendering(webgl, renderable);
        render(minDistTex, true);

        setupGroupIdRendering(webgl, renderable);
        render(texture, false);
    }

    // printTextureImage(readTexture(webgl, minDistTex), { scale: 0.75 });

    return { texture, scale, bbox: expandedBox, gridDim: dim, gridTexDim, gridTexScale, radiusFactor, resolution, maxRadius };
}

function calcGaussianDensityTexture3d(webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, texture?: Texture): _GaussianDensityTextureData {
    // console.log('3d');
    const { gl, resources, state, extensions: { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } } = webgl;
    const { smoothness, resolution } = props;

    const { drawCount, positions, radii, groups, scale, expandedBox, dim, maxRadius } = prepareGaussianDensityData(position, box, radius, props);
    const [dx, dy, dz] = dim;

    const minDistTex = getTexture('min-dist-3d', webgl, 'volume-uint8', 'rgba', 'ubyte', 'nearest');
    minDistTex.define(dx, dy, dz);

    const gridTexScale = Vec2.create(1, 1);
    const radiusFactor = maxRadius * 2;

    const renderable = getGaussianDensityRenderable(webgl, drawCount, positions, radii, groups, minDistTex, expandedBox, dim, dim, gridTexScale, smoothness, resolution, radiusFactor);

    //

    const { uCurrentSlice } = renderable.values;

    const framebuffer = getFramebuffer(webgl);
    framebuffer.bind();
    setRenderingDefaults(webgl);
    state.viewport(0, 0, dx, dy);
    state.scissor(0, 0, dx, dy);

    if (!texture) texture = colorBufferHalfFloat && textureHalfFloat
        ? resources.texture('volume-float16', 'rgba', 'fp16', 'linear')
        : colorBufferFloat && textureFloat
            ? resources.texture('volume-float32', 'rgba', 'float', 'linear')
            : resources.texture('volume-uint8', 'rgba', 'ubyte', 'linear');
    texture.define(dx, dy, dz);

    function render(fbTex: Texture, clear: boolean) {
        state.currentRenderItemId = -1;
        for (let i = 0; i < dz; ++i) {
            ValueCell.update(uCurrentSlice, i);
            fbTex.attachFramebuffer(framebuffer, 0, i);
            if (clear) gl.clear(gl.COLOR_BUFFER_BIT);
            renderable.render();
        }
        gl.flush();
    }

    setupDensityRendering(webgl, renderable);
    render(texture, true);

    setupMinDistanceRendering(webgl, renderable);
    render(minDistTex, true);

    setupGroupIdRendering(webgl, renderable);
    render(texture, false);

    return { texture, scale, bbox: expandedBox, gridDim: dim, gridTexDim: dim, gridTexScale, radiusFactor, resolution, maxRadius };
}

//

function prepareGaussianDensityData(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps) {
    const { resolution, radiusOffset } = props;
    const scaleFactor = 1 / resolution;

    const { indices, x, y, z, id } = position;
    const n = OrderedSet.size(indices);

    const positions = new Float32Array(n * 3);
    const radii = new Float32Array(n);
    const groups = new Float32Array(n);

    let maxRadius = 0;

    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i);

        positions[i * 3] = x[j];
        positions[i * 3 + 1] = y[j];
        positions[i * 3 + 2] = z[j];
        const r = radius(j) + radiusOffset;
        if (maxRadius < r) maxRadius = r;
        radii[i] = r;
        groups[i] = id ? id[i] : i;
    }

    const pad = maxRadius * 2 + resolution * 4;
    const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));
    const scaledBox = Box3D.scale(Box3D(), expandedBox, scaleFactor);
    const dim = Box3D.size(Vec3(), scaledBox);
    Vec3.ceil(dim, dim);

    const scale = Vec3.create(resolution, resolution, resolution);

    return { drawCount: n, positions, radii, groups, scale, expandedBox, dim, maxRadius };
}

function getGaussianDensityRenderable(webgl: WebGLContext, drawCount: number, positions: Float32Array, radii: Float32Array, groups: Float32Array, minDistanceTexture: Texture, box: Box3D, gridDim: Vec3, gridTexDim: Vec3, gridTexScale: Vec2, smoothness: number, resolution: number, radiusFactor: number) {
    // console.log('radiusFactor', radiusFactor);
    if (webgl.namedComputeRenderables[GaussianDensityName]) {
        const extent = Vec3.sub(Vec3(), box.max, box.min);
        const v = webgl.namedComputeRenderables[GaussianDensityName].values;

        ValueCell.updateIfChanged(v.drawCount, drawCount);
        ValueCell.updateIfChanged(v.instanceCount, 1);

        ValueCell.update(v.aRadius, radii);
        ValueCell.update(v.aPosition, positions);
        ValueCell.update(v.aGroup, groups);

        ValueCell.updateIfChanged(v.uCurrentSlice, 0);
        ValueCell.updateIfChanged(v.uCurrentX, 0);
        ValueCell.updateIfChanged(v.uCurrentY, 0);
        ValueCell.update(v.uBboxMin, box.min);
        ValueCell.update(v.uBboxSize, extent);
        ValueCell.update(v.uGridDim, gridDim);
        ValueCell.update(v.uGridTexDim, gridTexDim);
        ValueCell.update(v.uGridTexScale, gridTexScale);
        ValueCell.updateIfChanged(v.uAlpha, smoothness);
        ValueCell.updateIfChanged(v.uResolution, resolution);
        ValueCell.updateIfChanged(v.uRadiusFactorInv, 1 / radiusFactor);
        ValueCell.update(v.tMinDistanceTex, minDistanceTexture);

        ValueCell.updateIfChanged(v.dGridTexType, minDistanceTexture.getDepth() > 0 ? '3d' : '2d');
        ValueCell.updateIfChanged(v.dCalcType, 'density');

        webgl.namedComputeRenderables[GaussianDensityName].update();
    } else {
        webgl.namedComputeRenderables[GaussianDensityName] = createGaussianDensityRenderable(webgl, drawCount, positions, radii, groups, minDistanceTexture, box, gridDim, gridTexDim, gridTexScale, smoothness, resolution, radiusFactor);
    }
    return webgl.namedComputeRenderables[GaussianDensityName];
}

function createGaussianDensityRenderable(webgl: WebGLContext, drawCount: number, positions: Float32Array, radii: Float32Array, groups: Float32Array, minDistanceTexture: Texture, box: Box3D, gridDim: Vec3, gridTexDim: Vec3, gridTexScale: Vec2, smoothness: number, resolution: number, radiusFactor: number) {
    const extent = Vec3.sub(Vec3(), box.max, box.min);

    const values: GaussianDensityValues = {
        drawCount: ValueCell.create(drawCount),
        instanceCount: ValueCell.create(1),

        aRadius: ValueCell.create(radii),
        aPosition: ValueCell.create(positions),
        aGroup: ValueCell.create(groups),

        uCurrentSlice: ValueCell.create(0),
        uCurrentX: ValueCell.create(0),
        uCurrentY: ValueCell.create(0),
        uBboxMin: ValueCell.create(box.min),
        uBboxSize: ValueCell.create(extent),
        uGridDim: ValueCell.create(gridDim),
        uGridTexDim: ValueCell.create(gridTexDim),
        uGridTexScale: ValueCell.create(gridTexScale),
        uAlpha: ValueCell.create(smoothness),
        uResolution: ValueCell.create(resolution),
        uRadiusFactorInv: ValueCell.create(1 / radiusFactor),
        tMinDistanceTex: ValueCell.create(minDistanceTexture),

        dGridTexType: ValueCell.create(minDistanceTexture.getDepth() > 0 ? '3d' : '2d'),
        dCalcType: ValueCell.create('density'),
    };

    const schema = { ...GaussianDensitySchema };
    const shaderCode = ShaderCode(GaussianDensityName, gaussianDensity_vert, gaussianDensity_frag);
    const renderItem = createComputeRenderItem(webgl, 'points', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function setRenderingDefaults(ctx: WebGLContext) {
    const { gl, state } = ctx;
    state.disable(gl.CULL_FACE);
    state.enable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.enable(gl.SCISSOR_TEST);
    state.depthMask(false);
    state.clearColor(0, 0, 0, 0);
}

function setupMinDistanceRendering(webgl: WebGLContext, renderable: GaussianDensityRenderable) {
    const { gl, state } = webgl;
    ValueCell.update(renderable.values.dCalcType, 'minDistance');
    renderable.update();
    state.colorMask(false, false, false, true);
    state.blendFunc(gl.ONE, gl.ONE);
    // the shader writes 1 - dist so we set blending to MAX
    if (!webgl.extensions.blendMinMax) {
        throw new Error('GPU gaussian surface calculation requires EXT_blend_minmax');
    }
    state.blendEquation(webgl.extensions.blendMinMax.MAX);
}

function setupDensityRendering(webgl: WebGLContext, renderable: GaussianDensityRenderable) {
    const { gl, state } = webgl;
    ValueCell.update(renderable.values.dCalcType, 'density');
    renderable.update();
    state.colorMask(false, false, false, true);
    state.blendFunc(gl.ONE, gl.ONE);
    state.blendEquation(gl.FUNC_ADD);
}

function setupGroupIdRendering(webgl: WebGLContext, renderable: GaussianDensityRenderable) {
    const { gl, state } = webgl;
    ValueCell.update(renderable.values.dCalcType, 'groupId');
    renderable.update();
    // overwrite color, don't change alpha
    state.colorMask(true, true, true, false);
    state.blendFunc(gl.ONE, gl.ZERO);
    state.blendEquation(gl.FUNC_ADD);
}

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

function fieldFromTexture2d(ctx: WebGLContext, texture: Texture, dim: Vec3, texDim: Vec3) {
    // console.time('fieldFromTexture2d')
    const [dx, dy, dz] = dim;
    const [width, height] = texDim;
    const fboTexCols = Math.floor(width / dx);

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array);
    const data = space.create();
    const field = Tensor.create(space, data);
    const idData = space.create();
    const idField = Tensor.create(space, idData);

    const image = new Uint8Array(width * height * 4);

    const framebuffer = getFramebuffer(ctx);
    framebuffer.bind();
    texture.attachFramebuffer(framebuffer, 0);
    ctx.readPixels(0, 0, width, height, image);

    // printImageData(createImageData(image, width, height), 1/3)

    let j = 0;
    let tmpCol = 0;
    let tmpRow = 0;
    for (let iz = 0; iz < dz; ++iz) {
        if (tmpCol >= fboTexCols) {
            tmpCol = 0;
            tmpRow += dy;
        }
        for (let iy = 0; iy < dy; ++iy) {
            for (let ix = 0; ix < dx; ++ix) {
                const idx = 4 * (tmpCol * dx + (iy + tmpRow) * width + ix);
                data[j] = image[idx + 3] / 255;
                idData[j] = unpackRGBToInt(image[idx], image[idx + 1], image[idx + 2]);
                j++;
            }
        }
        tmpCol++;
    }

    // console.timeEnd('fieldFromTexture2d')

    return { field, idField };
}