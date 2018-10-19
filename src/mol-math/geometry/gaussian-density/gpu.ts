/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import { RuntimeContext } from 'mol-task'
import { PositionData, DensityData, DensityTextureData } from '../common'
import { Box3D } from '../../geometry'
import { GaussianDensityProps, getDelta } from '../gaussian-density'
import { OrderedSet } from 'mol-data/int'
import { Vec3, Tensor, Mat4 } from '../../linear-algebra'
import { GaussianDensityValues } from 'mol-gl/renderable/gaussian-density'
import { ValueCell, defaults } from 'mol-util'
import { RenderableState } from 'mol-gl/renderable'
import { createRenderable, createGaussianDensityRenderObject } from 'mol-gl/render-object'
import { Context, createContext, getGLContext } from 'mol-gl/webgl/context';
import { createFramebuffer } from 'mol-gl/webgl/framebuffer';
import { createTexture, Texture } from 'mol-gl/webgl/texture';
import { GLRenderingContext } from 'mol-gl/webgl/compat';
import { decodeIdRGBA } from 'mol-geo/geometry/picking';

export async function GaussianDensityGPU(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps): Promise<DensityData> {
    const webgl = defaults(props.webgl, getWebGLContext())

    const { transform, texture, gridDimension } = await GaussianDensityTexture(ctx, webgl, position, box, radius, props)

    const { field, idField } = webgl.isWebGL2 ?
        fieldFromTexture3d(webgl, texture, gridDimension) :
        fieldFromTexture2d(webgl, texture, gridDimension)

    console.log(idField)

    return { field, idField, transform }
}

export async function GaussianDensityTexture(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, oldTexture?: Texture): Promise<DensityTextureData> {
    console.time(`GaussianDensityTexture, ${webgl.isWebGL2 ? '3d' : '2d'}`)
    const { texture, scale, bbox, dim } = webgl.isWebGL2 ?
        await GaussianDensityTexture3d(ctx, webgl, position, box, radius, props, oldTexture) :
        await GaussianDensityTexture2d(ctx, webgl, position, box, radius, props, oldTexture)
    console.timeEnd(`GaussianDensityTexture, ${webgl.isWebGL2 ? '3d' : '2d'}`)

    const transform = Mat4.identity()
    Mat4.fromScaling(transform, scale)
    Mat4.setTranslation(transform, bbox.min)

    return { transform, texture, bbox, gridDimension: dim }
}

//

async function GaussianDensityTexture2d(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, texture?: Texture) {
    const { smoothness } = props

    const { drawCount, positions, radii, groups, delta, expandedBox, dim } = await prepareGaussianDensityData(ctx, position, box, radius, props)
    const [ dx, dy, dz ] = dim
    const minDistanceTexture = getMinDistanceTexture(webgl, dx, dy)

    const renderObject = getGaussianDensityRenderObject(webgl, drawCount, positions, radii, groups, minDistanceTexture, expandedBox, dim, smoothness)
    const renderable = createRenderable(webgl, renderObject)

    //

    const maxTexSize = webgl.maxTextureSize
    let fboTexDimX = 0
    let fboTexDimY = dim[1]
    let fboTexRows = 1
    let fboTexCols = dim[0]
    if (maxTexSize < dim[0] * dim[2]) {
        fboTexCols =  Math.floor(maxTexSize / dim[0])
        fboTexRows = Math.ceil(dim[2] / fboTexCols)
        fboTexDimX = fboTexCols * dim[0]
        fboTexDimY *= fboTexRows
    } else {
        fboTexDimX = dim[0] * dim[2]
    }

    //

    const { gl } = webgl
    const { uCurrentSlice, uCurrentX, uCurrentY } = renderObject.values

    const framebuffer = createFramebuffer(webgl)
    framebuffer.bind()

    if (!texture) {
        texture = createTexture(webgl, 'image-uint8', 'rgba', 'ubyte', 'linear')
    }
    texture.define(fboTexDimX, fboTexDimY)

    const program = renderable.getProgram('draw')
    program.use()
    setRenderingDefaults(gl)
    texture.attachFramebuffer(framebuffer, 0)

    let currCol = 0
    let currY = 0
    let currX = 0
    for (let i = 0; i < dz; ++i) {
        if (currCol >= fboTexCols) {
            currCol -= fboTexCols
            currY += dy
            currX = 0
        }
        gl.viewport(currX, currY, dx, dy)
        ValueCell.update(uCurrentSlice, i)
        ValueCell.update(uCurrentX, currX)
        ValueCell.update(uCurrentY, currY)
        renderable.render('draw')
        ++currCol
        currX += dx
    }

    framebuffer.destroy() // clean up

    await ctx.update({ message: 'gpu gaussian density calculation' });
    await webgl.waitForGpuCommandsComplete()

    return { texture, scale: Vec3.inverse(Vec3.zero(), delta), bbox: expandedBox, dim }
}

async function GaussianDensityTexture3d(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, texture?: Texture) {
    const { smoothness } = props

    const { drawCount, positions, radii, groups, delta, expandedBox, dim } = await prepareGaussianDensityData(ctx, position, box, radius, props)
    const [ dx, dy, dz ] = dim
    const minDistanceTexture = createTexture(webgl, 'volume-uint8', 'rgba', 'ubyte', 'nearest')
    minDistanceTexture.define(dx, dy, dz)

    const renderObject = getGaussianDensityRenderObject(webgl, drawCount, positions, radii, groups, minDistanceTexture, expandedBox, dim, smoothness)
    const renderable = createRenderable(webgl, renderObject)

    //

    const gl = webgl.gl as WebGL2RenderingContext
    const { uCurrentSlice } = renderObject.values

    const framebuffer = createFramebuffer(webgl)
    framebuffer.bind()

    gl.viewport(0, 0, dx, dy)
    setRenderingDefaults(gl)

    if (!texture) {
        texture = createTexture(webgl, 'volume-uint8', 'rgba', 'ubyte', 'linear')
    }
    texture.define(dx, dy, dz)

    ValueCell.update(renderable.values.dCalcType, 'minDistance')
    renderable.update()
    const programMinDistance = renderable.getProgram('draw')
    programMinDistance.use()
    gl.blendFunc(gl.ONE, gl.ONE)
    gl.blendEquation(gl.MAX)
    for (let i = 0; i < dz; ++i) {
        ValueCell.update(uCurrentSlice, i)
        minDistanceTexture.attachFramebuffer(framebuffer, 0, i)
        renderable.render('draw')
    }

    ValueCell.update(renderable.values.dCalcType, 'density')
    renderable.update()
    const programDensity = renderable.getProgram('draw')
    programDensity.use()
    gl.blendFunc(gl.ONE, gl.ONE)
    gl.blendEquation(gl.FUNC_ADD)
    for (let i = 0; i < dz; ++i) {
        ValueCell.update(uCurrentSlice, i)
        texture.attachFramebuffer(framebuffer, 0, i)
        renderable.render('draw')
    }

    ValueCell.update(renderable.values.dCalcType, 'groupId')
    renderable.update()
    const programGroupId = renderable.getProgram('draw')
    programGroupId.use()
    gl.blendFuncSeparate(gl.ONE, gl.ZERO, gl.ZERO, gl.ONE)
    gl.blendEquation(gl.FUNC_ADD)
    for (let i = 0; i < dz; ++i) {
        ValueCell.update(uCurrentSlice, i)
        texture.attachFramebuffer(framebuffer, 0, i)
        renderable.render('draw')
    }

    framebuffer.destroy() // clean up

    await ctx.update({ message: 'gpu gaussian density calculation' });
    await webgl.waitForGpuCommandsComplete()

    return { texture, scale: Vec3.inverse(Vec3.zero(), delta), bbox: expandedBox, dim }
}

//

let webglContext: Context
function getWebGLContext() {
    if (webglContext) return webglContext
    const canvas = document.createElement('canvas')
    const gl = getGLContext(canvas, {
        alpha: true,
        antialias: false,
        depth: false,
        preserveDrawingBuffer: true
    })
    if (!gl) throw new Error('Could not create a WebGL rendering context')
    webglContext = createContext(gl)
    return webglContext
}

async function prepareGaussianDensityData(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps) {
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
    console.log('grid dim gpu', dim)

    return { drawCount: n, positions, radii, groups, delta, expandedBox, dim }
}

function getMinDistanceTexture(webgl: Context, width: number, height: number) {
    const minDistanceTexture = createTexture(webgl, 'image-uint8', 'rgba', 'ubyte', 'nearest')
    minDistanceTexture.define(width, height)
    return minDistanceTexture
}

function getGaussianDensityRenderObject(webgl: Context, drawCount: number, positions: Float32Array, radii: Float32Array, groups: Float32Array, minDistanceTexture: Texture, box: Box3D, dimensions: Vec3, smoothness: number) {
    const extent = Vec3.sub(Vec3.zero(), box.max, box.min)

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
        uBboxMax: ValueCell.create(box.max),
        uBboxSize: ValueCell.create(extent),
        uGridDim: ValueCell.create(dimensions),
        uAlpha: ValueCell.create(smoothness),
        tMinDistanceTex: ValueCell.create(minDistanceTexture),

        dCalcType: ValueCell.create('density'),
    }
    const state: RenderableState = {
        visible: true,
        depthMask: false
    }

    const renderObject = createGaussianDensityRenderObject(values, state)

    return renderObject
}

function setRenderingDefaults(gl: GLRenderingContext) {
    gl.disable(gl.CULL_FACE)
    gl.frontFace(gl.CCW)
    gl.cullFace(gl.BACK)

    gl.blendFunc(gl.ONE, gl.ONE)
    gl.blendEquation(gl.FUNC_ADD)
    gl.enable(gl.BLEND)
}

function fieldFromTexture2d(ctx: Context, texture: Texture, dim: Vec3) {
    console.time('fieldFromTexture2d')
    const { gl } = ctx
    const [ dx, dy, dz ] = dim
    const { width, height } = texture
    const fboTexCols = Math.floor(width / dx)

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)
    const idData = space.create()
    const idField = Tensor.create(space, idData)

    const image = new Uint8Array(width * height * 4)

    const framebuffer = createFramebuffer(ctx)
    framebuffer.bind()

    texture.attachFramebuffer(framebuffer, 0)
    gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, image)

    let idx = 0
    let tmpCol = 0
    let tmpRow = 0
    for (let iz = 0; iz < dz; ++iz) {
        if (tmpCol >= fboTexCols ) {
            tmpCol = 0
            tmpRow += dy
        }
        for (let iy = 0; iy < dy; ++iy) {
            for (let ix = 0; ix < dx; ++ix) {
                data[idx] = image[4 * (tmpCol * dx + (iy + tmpRow) * width + ix) + 3] / 255
                idx++
            }
        }
        tmpCol++
    }

    framebuffer.destroy()
    console.timeEnd('fieldFromTexture2d')

    return { field, idField }
}

function fieldFromTexture3d(ctx: Context, texture: Texture, dim: Vec3) {
    console.time('fieldFromTexture3d')
    const { gl } = ctx
    const { width, height, depth } = texture

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)
    const idData = space.create()
    const idField = Tensor.create(space, idData)

    const slice = new Uint8Array(width * height * 4)

    const framebuffer = createFramebuffer(ctx)
    framebuffer.bind()

    let j = 0
    for (let i = 0; i < depth; ++i) {
        texture.attachFramebuffer(framebuffer, 0, i)
        gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, slice)
        for (let iy = 0; iy < height; ++iy) {
            for (let ix = 0; ix < width; ++ix) {
                const idx = 4 * (iy * width + ix)
                data[j] = slice[idx + 3] / 255
                idData[j] = decodeIdRGBA(slice[idx], slice[idx + 1], slice[idx + 2])
                ++j
            }
        }
    }

    framebuffer.destroy()
    console.timeEnd('fieldFromTexture3d')

    return { field, idField }
}