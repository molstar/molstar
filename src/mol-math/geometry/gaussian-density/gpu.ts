/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import { RuntimeContext } from 'mol-task'
import { PositionData, DensityData } from '../common'
import { Box3D } from '../../geometry'
import { GaussianDensityProps, getDelta } from '../gaussian-density'
import { OrderedSet } from 'mol-data/int'
import { Vec3, Tensor, Mat4 } from '../../linear-algebra'
import { GaussianDensityValues } from 'mol-gl/renderable/gaussian-density'
import { ValueCell } from 'mol-util'
import { RenderableState } from 'mol-gl/renderable'
import { createRenderable, createGaussianDensityRenderObject } from 'mol-gl/render-object'
import { createRenderTarget } from 'mol-gl/webgl/render-target'
import { Context, createContext, getGLContext } from 'mol-gl/webgl/context';
import { createFramebuffer } from 'mol-gl/webgl/framebuffer';
import { createTexture, Texture, TextureAttachment } from 'mol-gl/webgl/texture';
import { GLRenderingContext } from 'mol-gl/webgl/compat';

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

export async function GaussianDensityGPU(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps): Promise<DensityData> {
    // TODO allow passing a context via props
    const webgl = getWebGLContext()

    if (webgl.maxDrawBuffers > 0) {
        console.log('GaussianDensityMultiDrawBuffer')
        const { texture, scale, bbox, dim } = await GaussianDensityMultiDrawBuffer(ctx, webgl, position, box, radius, props)

        console.time('gpu gaussian density 3d texture read')
        const field = fieldFromTexture3d(webgl, texture, dim)
        console.timeEnd('gpu gaussian density 3d texture read')

        const idData = field.space.create()
        const idField = Tensor.create(field.space, idData)

        const transform = Mat4.identity()
        Mat4.fromScaling(transform, scale)
        Mat4.setTranslation(transform, bbox.min)

        const renderTarget = createRenderTarget(webgl, dim[0], dim[1])

        return { field, idField, transform, renderTarget, bbox, gridDimension: dim }
    } else {
        console.log('GaussianDensitySingleDrawBuffer')
        return GaussianDensitySingleDrawBuffer(ctx, webgl, position, box, radius, props)
    }
}

async function prepareGaussianDensityData(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps) {
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

//

function setRenderingDefaults(gl: GLRenderingContext) {
    gl.disable(gl.CULL_FACE)
    gl.frontFace(gl.CCW)
    gl.cullFace(gl.BACK)

    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
    gl.blendEquation(gl.FUNC_ADD)
    gl.enable(gl.BLEND)
}

//

async function GaussianDensitySingleDrawBuffer(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps): Promise<DensityData> {
    const { readSlices, smoothness } = props

    const { drawCount, positions, radii, delta, expandedBox, dim } = await prepareGaussianDensityData(ctx, position, box, radius, props)
    const renderObject = getGaussianDensityRenderObject(webgl, drawCount, positions, radii, expandedBox, dim, smoothness)
    const renderable = createRenderable(webgl, renderObject)

    //

    // TODO fallback to lower resolution when texture size is not large enough
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

    console.log('dim', dim, 'cols', fboTexCols, 'rows', fboTexRows)

    //

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    const idData = space.create()
    const idField = Tensor.create(space, idData)

    //

    const { gl } = webgl
    const { uCurrentSlice, uCurrentX, uCurrentY } = renderObject.values

    const program = renderable.getProgram('draw')
    const renderTarget = createRenderTarget(webgl, fboTexDimX, fboTexDimY)

    program.use()
    renderTarget.bind()
    setRenderingDefaults(gl)

    const slice = new Uint8Array(dim[0] * dim[1] * 4)

    console.time('gpu gaussian density slices')
    let currCol = 0
    let currY = 0
    let currX = 0
    let j = 0
    for (let i = 0; i < dim[2]; ++i) {
        if (currCol >= fboTexCols) {
            currCol -= fboTexCols
            currY += dim[1]
            currX = 0
        }
        gl.viewport(currX, currY, dim[0], dim[1])
        ValueCell.update(uCurrentSlice, i)
        ValueCell.update(uCurrentX, currX)
        ValueCell.update(uCurrentY, currY)
        renderable.render('draw')
        if (readSlices) {
            renderTarget.readBuffer(currX, currY, dim[0], dim[1], slice)
            for (let iy = 0; iy < dim[1]; ++iy) {
                for (let ix = 0; ix < dim[0]; ++ix) {
                    data[j] = slice[4 * (iy * dim[0] + ix)] / 255
                    ++j
                }
            }
        }
        ++currCol
        currX += dim[0]
    }
    console.timeEnd('gpu gaussian density slices')

    //

    if (!readSlices) {
        console.time('gpu gaussian density full')
        renderTarget.getBuffer()
        const { array } = renderTarget.image
        let idx = 0
        let tmpCol = 0
        let tmpRow = 0
        for (let iz = 0; iz < dim[2]; ++iz) {
            if (tmpCol >= fboTexCols ) {
                tmpCol = 0
                tmpRow += dim[1]
            }
            for (let iy = 0; iy < dim[1]; ++iy) {
                for (let ix = 0; ix < dim[0]; ++ix) {
                    data[idx] = array[4 * (tmpCol * dim[0] + (iy + tmpRow) * fboTexDimX + ix)] / 255
                    idx++
                }
            }
            tmpCol++
        }
        console.timeEnd('gpu gaussian density full')
    }

    //

    const transform = Mat4.identity()
    Mat4.fromScaling(transform, Vec3.inverse(Vec3.zero(), delta))
    Mat4.setTranslation(transform, expandedBox.min)

    return { field, idField, transform, renderTarget, bbox: expandedBox, gridDimension: dim }
}

async function GaussianDensityMultiDrawBuffer(ctx: RuntimeContext, webgl: Context, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps) {
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

    const texture = createTexture(webgl, 'volume-uint8', 'rgba', 'ubyte', 'linear')
    texture.define(dx, dy, dz)

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

    gl.viewport(0, 0, dx, dy)
    setRenderingDefaults(gl)

    //

    const dzMulti = Math.floor(dz / drawBuffers) * drawBuffers

    const programMulti = renderable.getProgram('draw')
    programMulti.use()

    console.time('gpu gaussian density 3d texture multi')
    for (let i = 0; i < dzMulti; i += drawBuffers) {
        ValueCell.update(uCurrentSlice, i)
        for (let k = 0; k < drawBuffers; ++k) {
            texture.attachFramebuffer(framebuffer, k as TextureAttachment, i + k)
        }
        renderable.render('draw')
    }
    console.timeEnd('gpu gaussian density 3d texture multi')

    ValueCell.updateIfChanged(renderable.values.dDrawBuffers, 1)
    renderable.update()
    const programSingle = renderable.getProgram('draw')
    programSingle.use()

    console.time('gpu gaussian density 3d texture single')
    for (let i = dzMulti; i < dz; ++i) {
        ValueCell.update(uCurrentSlice, i)
        texture.attachFramebuffer(framebuffer, 0, i)
        renderable.render('draw')
    }
    console.timeEnd('gpu gaussian density 3d texture single')

    // must detach framebuffer attachments before reading is possible
    for (let k = 0; k < drawBuffers; ++k) {
        texture.detachFramebuffer(framebuffer, k as TextureAttachment)
    }
    framebuffer.destroy() // clean up

    // throw new Error('foo')

    return { texture, scale: Vec3.inverse(Vec3.zero(), delta), bbox: expandedBox, dim }
}

//

function fieldFromTexture3d(ctx: Context, texture: Texture, dim: Vec3) {
    const { gl } = ctx
    const [ dx, dy, dz ] = dim

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    const slice = new Uint8Array(dx * dy * 4)

    const framebuffer = createFramebuffer(ctx)
    framebuffer.bind()

    let j = 0
    for (let i = 0; i < dz; ++i) {
        texture.attachFramebuffer(framebuffer, 0, i)
        gl.readPixels(0, 0, dx, dy, gl.RGBA, gl.UNSIGNED_BYTE, slice)
        for (let iy = 0; iy < dim[1]; ++iy) {
            for (let ix = 0; ix < dim[0]; ++ix) {
                data[j] = slice[4 * (iy * dim[0] + ix)] / 255
                ++j
            }
        }
    }

    framebuffer.destroy()

    return field
}