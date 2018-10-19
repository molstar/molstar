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
import { createTexture, Texture, TextureAttachment } from 'mol-gl/webgl/texture';
import { GLRenderingContext } from 'mol-gl/webgl/compat';

export async function GaussianDensityGPU(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps): Promise<DensityData> {
    const webgl = defaults(props.webgl, getWebGLContext())

    const { transform, texture, gridDimension } = await GaussianDensityTexture(ctx, webgl, position, box, radius, props)

    const field = webgl.maxDrawBuffers > 0 ?
        fieldFromTexture3d(webgl, texture, gridDimension) :
        fieldFromTexture2d(webgl, texture, gridDimension)

    const idData = field.space.create()
    const idField = Tensor.create(field.space, idData)

    return { field, idField, transform }
}

export async function GaussianDensityTexture(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, oldTexture?: Texture): Promise<DensityTextureData> {
    console.time(`GaussianDensityTexture, ${webgl.maxDrawBuffers > 0 ? 'multi' : 'single'}`)
    const { texture, scale, bbox, dim } = webgl.maxDrawBuffers > 0 ?
        await GaussianDensityMultiDrawBuffer(ctx, webgl, position, box, radius, props, oldTexture) :
        await GaussianDensitySingleDrawBuffer(ctx, webgl, position, box, radius, props, oldTexture)
    console.timeEnd(`GaussianDensityTexture, ${webgl.maxDrawBuffers > 0 ? 'multi' : 'single'}`)

    const transform = Mat4.identity()
    Mat4.fromScaling(transform, scale)
    Mat4.setTranslation(transform, bbox.min)

    return { transform, texture, bbox, gridDimension: dim }
}

//

async function GaussianDensitySingleDrawBuffer(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, texture?: Texture) {
    const { smoothness } = props

    const { drawCount, positions, radii, delta, expandedBox, dim } = await prepareGaussianDensityData(ctx, position, box, radius, props)
    const [ dx, dy, dz ] = dim
    const renderObject = getGaussianDensityRenderObject(webgl, drawCount, positions, radii, expandedBox, dim, smoothness)
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

async function GaussianDensityMultiDrawBuffer(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, texture?: Texture) {
    const { smoothness } = props

    const { drawCount, positions, radii, delta, expandedBox, dim } = await prepareGaussianDensityData(ctx, position, box, radius, props)
    const [ dx, dy, dz ] = dim
    const renderObject = getGaussianDensityRenderObject(webgl, drawCount, positions, radii, expandedBox, dim, smoothness)
    const renderable = createRenderable(webgl, renderObject)
    const drawBuffers = Math.min(8, webgl.maxDrawBuffers)

    //

    const gl = webgl.gl as WebGL2RenderingContext
    const { uCurrentSlice } = renderObject.values

    const framebuffer = createFramebuffer(webgl)
    framebuffer.bind()

    setDrawBuffers(gl, drawBuffers)
    gl.viewport(0, 0, dx, dy)
    setRenderingDefaults(gl)

    if (!texture) {
        texture = createTexture(webgl, 'volume-uint8', 'rgba', 'ubyte', 'linear')
    }
    texture.define(dx, dy, dz)

    // z-slices to be render with multi render targets
    const dzMulti = Math.floor(dz / drawBuffers) * drawBuffers

    // render multi target
    const programMulti = renderable.getProgram('draw')
    programMulti.use()
    for (let i = 0; i < dzMulti; i += drawBuffers) {
        ValueCell.update(uCurrentSlice, i)
        for (let k = 0; k < drawBuffers; ++k) {
            texture.attachFramebuffer(framebuffer, k as TextureAttachment, i + k)
        }
        renderable.render('draw');
    }

    // render single target
    ValueCell.updateIfChanged(renderable.values.dDrawBuffers, 1)
    renderable.update()
    const programSingle = renderable.getProgram('draw')
    programSingle.use()
    for (let i = dzMulti; i < dz; ++i) {
        ValueCell.update(uCurrentSlice, i)
        texture.attachFramebuffer(framebuffer, 0, i)
        renderable.render('draw')
    }

    // must detach framebuffer attachments before reading is possible
    for (let k = 0; k < drawBuffers; ++k) {
        texture.detachFramebuffer(framebuffer, k as TextureAttachment)
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

    let maxRadius = 0

    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i);

        positions[i * 3] = x[j]
        positions[i * 3 + 1] = y[j]
        positions[i * 3 + 2] = z[j]
        const r = radius(j) + radiusOffset
        if (maxRadius < r) maxRadius = r
        radii[i] = r

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

    return { drawCount: n, positions, radii, delta, expandedBox, dim }
}

function getGaussianDensityRenderObject(webgl: Context, drawCount: number, positions: Float32Array, radii: Float32Array, box: Box3D, dimensions: Vec3, smoothness: number) {
    const extent = Vec3.sub(Vec3.zero(), box.max, box.min)

    const values: GaussianDensityValues = {
        drawCount: ValueCell.create(drawCount),
        instanceCount: ValueCell.create(1),

        aRadius: ValueCell.create(radii),
        aPosition: ValueCell.create(positions),

        uCurrentSlice: ValueCell.create(0),
        uCurrentX: ValueCell.create(0),
        uCurrentY: ValueCell.create(0),
        uBboxMin: ValueCell.create(box.min),
        uBboxMax: ValueCell.create(box.max),
        uBboxSize: ValueCell.create(extent),
        uGridDim: ValueCell.create(dimensions),
        uAlpha: ValueCell.create(smoothness),

        dDrawBuffers: ValueCell.create(Math.min(8, webgl.maxDrawBuffers)),
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

function setDrawBuffers(gl: WebGL2RenderingContext, drawBuffers: number) {
    if (drawBuffers === 1) {
        gl.drawBuffers([
            gl.COLOR_ATTACHMENT0,
        ]);
    } else if (drawBuffers === 4) {
        gl.drawBuffers([
            gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, gl.COLOR_ATTACHMENT2, gl.COLOR_ATTACHMENT3,
        ]);
    } else if (drawBuffers === 8) {
        gl.drawBuffers([
            gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, gl.COLOR_ATTACHMENT2, gl.COLOR_ATTACHMENT3,
            gl.COLOR_ATTACHMENT4, gl.COLOR_ATTACHMENT5, gl.COLOR_ATTACHMENT6, gl.COLOR_ATTACHMENT7,
        ]);
    }
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

    return field
}

function fieldFromTexture3d(ctx: Context, texture: Texture, dim: Vec3) {
    console.time('fieldFromTexture3d')
    const { gl } = ctx
    const { width, height, depth } = texture

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    const slice = new Uint8Array(width * height * 4)

    const framebuffer = createFramebuffer(ctx)
    framebuffer.bind()

    let j = 0
    for (let i = 0; i < depth; ++i) {
        texture.attachFramebuffer(framebuffer, 0, i)
        gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, slice)
        for (let iy = 0; iy < height; ++iy) {
            for (let ix = 0; ix < width; ++ix) {
                data[j] = slice[4 * (iy * width + ix) + 3] / 255
                ++j
            }
        }
    }

    framebuffer.destroy()
    console.timeEnd('fieldFromTexture3d')

    return field
}