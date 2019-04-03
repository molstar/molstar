/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import { RuntimeContext } from 'mol-task'
import { PositionData, DensityData, DensityTextureData } from '../common'
import { Box3D } from '../../geometry'
import { GaussianDensityGPUProps, getDelta } from '../gaussian-density'
import { OrderedSet } from 'mol-data/int'
import { Vec3, Tensor, Mat4 } from '../../linear-algebra'
import { ValueCell } from 'mol-util'
import { createComputeRenderable, ComputeRenderable } from 'mol-gl/renderable'
import { WebGLContext } from 'mol-gl/webgl/context';
import { createTexture, Texture } from 'mol-gl/webgl/texture';
import { GLRenderingContext } from 'mol-gl/webgl/compat';
import { decodeFloatRGB } from 'mol-util/float-packing';
import { ShaderCode } from 'mol-gl/shader-code';
import { createComputeRenderItem } from 'mol-gl/webgl/render-item';
import { ValueSpec, AttributeSpec, UniformSpec, TextureSpec, DefineSpec, Values } from 'mol-gl/renderable/schema';

export const GaussianDensitySchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    aRadius: AttributeSpec('float32', 1, 0),
    aPosition: AttributeSpec('float32', 3, 0),
    aGroup: AttributeSpec('float32', 1, 0),

    uCurrentSlice: UniformSpec('f'),
    uCurrentX: UniformSpec('f'),
    uCurrentY: UniformSpec('f'),
    uBboxMin: UniformSpec('v3'),
    uBboxMax: UniformSpec('v3'),
    uBboxSize: UniformSpec('v3'),
    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),
    uAlpha: UniformSpec('f'),
    tMinDistanceTex: TextureSpec('texture', 'rgba', 'float', 'nearest'),

    dGridTexType: DefineSpec('string', ['2d', '3d']),
    dCalcType: DefineSpec('string', ['density', 'minDistance', 'groupId']),
}

export const GaussianDensityShaderCode = ShaderCode(
    require('mol-gl/shader/gaussian-density.vert').default,
    require('mol-gl/shader/gaussian-density.frag').default,
    { standardDerivatives: false, fragDepth: false }
)

/** name for shared framebuffer used for gpu gaussian surface operations */
const FramebufferName = 'gaussian-density'

export async function GaussianDensityGPU(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, webgl: WebGLContext): Promise<DensityData> {
    // always use texture2d when the gaussian density needs to be downloaded from the GPU,
    // it's faster than texture3d
    // console.time('GaussianDensityTexture2d')
    const { scale, bbox, texture, dim } = await calcGaussianDensityTexture2d(ctx, webgl, position, box, radius, props)
    // webgl.waitForGpuCommandsCompleteSync()
    // console.timeEnd('GaussianDensityTexture2d')
    const { field, idField } = await fieldFromTexture2d(webgl, texture, dim)

    return { field, idField, transform: getTransform(scale, bbox) }
}

export async function GaussianDensityTexture(ctx: RuntimeContext, webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, oldTexture?: Texture): Promise<DensityTextureData> {
    return webgl.isWebGL2 ?
        await GaussianDensityTexture3d(ctx, webgl, position, box, radius, props, oldTexture) :
        await GaussianDensityTexture2d(ctx, webgl, position, box, radius, props, oldTexture)
}

export async function GaussianDensityTexture2d(ctx: RuntimeContext, webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, oldTexture?: Texture): Promise<DensityTextureData> {
    return finalizeGaussianDensityTexture(await calcGaussianDensityTexture2d(ctx, webgl, position, box, radius, props, oldTexture))
}

export async function GaussianDensityTexture3d(ctx: RuntimeContext, webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, oldTexture?: Texture): Promise<DensityTextureData> {
    return finalizeGaussianDensityTexture(await calcGaussianDensityTexture3d(ctx, webgl, position, box, radius, props, oldTexture))
}

function finalizeGaussianDensityTexture({ texture, scale, bbox, dim }: GaussianDensityTextureData): DensityTextureData {
    return { transform: getTransform(scale, bbox), texture, bbox, gridDimension: dim }
}

function getTransform(scale: Vec3, bbox: Box3D) {
    const transform = Mat4.identity()
    Mat4.fromScaling(transform, scale)
    Mat4.setTranslation(transform, bbox.min)
    return transform
}

//

type GaussianDensityTextureData = { texture: Texture, scale: Vec3, bbox: Box3D, dim: Vec3}

async function calcGaussianDensityTexture2d(ctx: RuntimeContext, webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, texture?: Texture): Promise<GaussianDensityTextureData> {
    const { smoothness } = props

    const { drawCount, positions, radii, groups, delta, expandedBox, dim } = await prepareGaussianDensityData(ctx, position, box, radius, props)
    const [ dx, dy, dz ] = dim
    const { texDimX, texDimY, texCols } = getTexture2dSize(dim)

    const minDistanceTexture = createTexture(webgl, 'image-float32', 'rgba', 'float', 'nearest')
    minDistanceTexture.define(texDimX, texDimY)

    const renderable = getGaussianDensityRenderable(webgl, drawCount, positions, radii, groups, minDistanceTexture, expandedBox, dim, smoothness)

    //

    const { gl, framebufferCache } = webgl
    const { uCurrentSlice, uCurrentX, uCurrentY } = renderable.values

    const framebuffer = framebufferCache.get(FramebufferName).value
    framebuffer.bind()
    setRenderingDefaults(gl)

    if (!texture) texture = createTexture(webgl, 'image-uint8', 'rgba', 'ubyte', 'linear')
    texture.define(texDimX, texDimY)

    function render(fbTex: Texture) {
        fbTex.attachFramebuffer(framebuffer, 0)
        let currCol = 0
        let currY = 0
        let currX = 0
        for (let i = 0; i < dz; ++i) {
            if (currCol >= texCols) {
                currCol -= texCols
                currY += dy
                currX = 0
            }
            gl.viewport(currX, currY, dx, dy)
            ValueCell.update(uCurrentSlice, i)
            ValueCell.update(uCurrentX, currX)
            ValueCell.update(uCurrentY, currY)
            renderable.render()
            ++currCol
            currX += dx
        }
    }

    setupMinDistanceRendering(webgl, renderable)
    render(minDistanceTexture)

    setupDensityRendering(webgl, renderable)
    render(texture)

    setupGroupIdRendering(webgl, renderable)
    render(texture)

    if (ctx.shouldUpdate) await ctx.update({ message: 'gpu gaussian density calculation' })

    return { texture, scale: Vec3.inverse(Vec3.zero(), delta), bbox: expandedBox, dim }
}

async function calcGaussianDensityTexture3d(ctx: RuntimeContext, webgl: WebGLContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, texture?: Texture): Promise<GaussianDensityTextureData> {
    const { smoothness } = props

    const { drawCount, positions, radii, groups, delta, expandedBox, dim } = await prepareGaussianDensityData(ctx, position, box, radius, props)
    const [ dx, dy, dz ] = dim
    const minDistanceTexture = createTexture(webgl, 'volume-float32', 'rgba', 'float', 'nearest')
    minDistanceTexture.define(dx, dy, dz)

    const renderable = getGaussianDensityRenderable(webgl, drawCount, positions, radii, groups, minDistanceTexture, expandedBox, dim, smoothness)

    //

    const { gl, framebufferCache } = webgl
    const { uCurrentSlice } = renderable.values

    const framebuffer = framebufferCache.get(FramebufferName).value
    framebuffer.bind()
    setRenderingDefaults(gl)
    gl.viewport(0, 0, dx, dy)

    if (!texture) texture = createTexture(webgl, 'volume-uint8', 'rgba', 'ubyte', 'linear')
    texture.define(dx, dy, dz)

    function render(fbTex: Texture) {
        for (let i = 0; i < dz; ++i) {
            ValueCell.update(uCurrentSlice, i)
            fbTex.attachFramebuffer(framebuffer, 0, i)
            renderable.render()
        }
    }

    setupMinDistanceRendering(webgl, renderable)
    render(minDistanceTexture)

    setupDensityRendering(webgl, renderable)
    render(texture)

    setupGroupIdRendering(webgl, renderable)
    render(texture)

    if (ctx.shouldUpdate) await ctx.update({ message: 'gpu gaussian density calculation' });

    return { texture, scale: Vec3.inverse(Vec3.zero(), delta), bbox: expandedBox, dim }
}

//

async function prepareGaussianDensityData(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps) {
    const { resolution, radiusOffset } = props

    const { indices, x, y, z } = position
    const n = OrderedSet.size(indices)

    const positions = new Float32Array(n * 3)
    const radii = new Float32Array(n)
    const groups = new Float32Array(n)

    let maxRadius = 0

    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i);

        positions[i * 3] = x[j]
        positions[i * 3 + 1] = y[j]
        positions[i * 3 + 2] = z[j]
        const r = radius(j) + radiusOffset
        if (maxRadius < r) maxRadius = r
        radii[i] = r
        groups[i] = i

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'preparing density data', current: i, max: n })
        }
    }

    const pad = maxRadius * 2 + resolution
    const expandedBox = Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad));
    const extent = Vec3.sub(Vec3.zero(), expandedBox.max, expandedBox.min)

    const delta = getDelta(expandedBox, resolution)
    const dim = Vec3.zero()
    Vec3.ceil(dim, Vec3.mul(dim, extent, delta))
    // console.log('grid dim gpu', dim)

    return { drawCount: n, positions, radii, groups, delta, expandedBox, dim }
}

function getGaussianDensityRenderable(webgl: WebGLContext, drawCount: number, positions: Float32Array, radii: Float32Array, groups: Float32Array, minDistanceTexture: Texture, box: Box3D, dimensions: Vec3, smoothness: number) {
    const extent = Vec3.sub(Vec3.zero(), box.max, box.min)
    const { texDimX, texDimY } = getTexture2dSize(dimensions)

    const values: Values<typeof GaussianDensitySchema> = {
        drawCount: ValueCell.create(drawCount),
        instanceCount: ValueCell.create(1),

        aRadius: ValueCell.create(radii),
        aPosition: ValueCell.create(positions),
        aGroup: ValueCell.create(groups),

        uCurrentSlice: ValueCell.create(0),
        uCurrentX: ValueCell.create(0),
        uCurrentY: ValueCell.create(0),
        uBboxMin: ValueCell.create(box.min),
        uBboxMax: ValueCell.create(box.max),
        uBboxSize: ValueCell.create(extent),
        uGridDim: ValueCell.create(dimensions),
        uGridTexDim: ValueCell.create(Vec3.create(texDimX, texDimY, 0)),
        uAlpha: ValueCell.create(smoothness),
        tMinDistanceTex: ValueCell.create(minDistanceTexture),

        dGridTexType: ValueCell.create(minDistanceTexture.depth > 0 ? '3d' : '2d'),
        dCalcType: ValueCell.create('minDistance'),
    }

    const schema = { ...GaussianDensitySchema }
    const shaderCode = GaussianDensityShaderCode
    const renderItem =  createComputeRenderItem(webgl, 'points', shaderCode, schema, values)

    return createComputeRenderable(renderItem, values)
}

function setRenderingDefaults(gl: GLRenderingContext) {
    gl.disable(gl.CULL_FACE)
    gl.frontFace(gl.CCW)
    gl.cullFace(gl.BACK)
    gl.enable(gl.BLEND)
}

function setupMinDistanceRendering(webgl: WebGLContext, renderable: ComputeRenderable<any>) {
    const { gl } = webgl
    ValueCell.update(renderable.values.dCalcType, 'minDistance')
    renderable.update()
    renderable.use()
    gl.blendFunc(gl.ONE, gl.ONE)
    // the shader writes 1 - dist so we set blending to MAX
    gl.blendEquation(webgl.extensions.blendMinMax.MAX)
}

function setupDensityRendering(webgl: WebGLContext, renderable: ComputeRenderable<any>) {
    const { gl } = webgl
    ValueCell.update(renderable.values.dCalcType, 'density')
    renderable.update()
    renderable.use()
    gl.blendFunc(gl.ONE, gl.ONE)
    gl.blendEquation(gl.FUNC_ADD)
}

function setupGroupIdRendering(webgl: WebGLContext, renderable: ComputeRenderable<any>) {
    const { gl } = webgl
    ValueCell.update(renderable.values.dCalcType, 'groupId')
    renderable.update()
    renderable.use()
    // overwrite color, don't change alpha
    gl.blendFuncSeparate(gl.ONE, gl.ZERO, gl.ZERO, gl.ONE)
    gl.blendEquation(gl.FUNC_ADD)
}

function getTexture2dSize(gridDim: Vec3) {
    const area = gridDim[0] * gridDim[1] * gridDim[2]
    const squareDim = Math.sqrt(area)
    const squareTexSize = Math.pow(2, Math.ceil(Math.log(squareDim) / Math.log(2)))

    let texDimX = 0
    let texDimY = gridDim[1]
    let texRows = 1
    let texCols = gridDim[2]
    if (squareTexSize < gridDim[0] * gridDim[2]) {
        texCols = Math.floor(squareTexSize / gridDim[0])
        texRows = Math.ceil(gridDim[2] / texCols)
        texDimX = texCols * gridDim[0]
        texDimY *= texRows
    } else {
        texDimX = gridDim[0] * gridDim[2]
    }
    return { texDimX, texDimY, texRows, texCols }
}

async function fieldFromTexture2d(ctx: WebGLContext, texture: Texture, dim: Vec3) {
    // console.time('fieldFromTexture2d')
    const { framebufferCache } = ctx
    const [ dx, dy, dz ] = dim
    const { width, height } = texture
    const fboTexCols = Math.floor(width / dx)

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)
    const idData = space.create()
    const idField = Tensor.create(space, idData)

    const image = new Uint8Array(width * height * 4)

    const framebuffer = framebufferCache.get(FramebufferName).value
    framebuffer.bind()
    texture.attachFramebuffer(framebuffer, 0)
    // TODO too slow, why? Too many checks if gpu ready???
    // await ctx.readPixelsAsync(0, 0, width, height, image)
    ctx.readPixels(0, 0, width, height, image)

    // printImageData(createImageData(image, width, height), 1/3)

    let j = 0
    let tmpCol = 0
    let tmpRow = 0
    for (let iz = 0; iz < dz; ++iz) {
        if (tmpCol >= fboTexCols ) {
            tmpCol = 0
            tmpRow += dy
        }
        for (let iy = 0; iy < dy; ++iy) {
            for (let ix = 0; ix < dx; ++ix) {
                const idx = 4 * (tmpCol * dx + (iy + tmpRow) * width + ix)
                data[j] = image[idx + 3] / 255
                idData[j] = decodeFloatRGB(image[idx], image[idx + 1], image[idx + 2])
                j++
            }
        }
        tmpCol++
    }

    // console.timeEnd('fieldFromTexture2d')

    return { field, idField }
}