/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '../geometry';
import { Vec3, Mat4, Tensor } from '../linear-algebra';
import { RuntimeContext, Task } from 'mol-task';
import { PositionData, DensityData } from './common';
import { OrderedSet } from 'mol-data/int';
import { createRenderable, createGaussianDensityRenderObject } from 'mol-gl/render-object';
import { createContext, Context } from 'mol-gl/webgl/context';
import { GaussianDensityValues } from 'mol-gl/renderable/gaussian-density';
import { RenderableState } from 'mol-gl/renderable';
import { ValueCell } from 'mol-util';
import { createRenderTarget } from 'mol-gl/webgl/render-target';

export const DefaultGaussianDensityProps = {
    resolution: 1,
    radiusOffset: 0,
    smoothness: 1.5,
    readSlices: false,
    useGpu: true,
}
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

function getDelta(box: Box3D, resolution: number) {
    const extent = Vec3.sub(Vec3.zero(), box.max, box.min)
    const size = Vec3.zero()
    Vec3.ceil(size, Vec3.scale(size, extent, resolution))
    const delta = Vec3.div(Vec3.zero(), extent, size)
    return delta
}

export function computeGaussianDensity(position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps) {
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, box, radius, props)
    });
}

export async function GaussianDensity(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps): Promise<DensityData> {
    if (props.useGpu) {
        return await GaussianDensityGPU(ctx, position, box, radius, props)
    } else {
        return await GaussianDensityCPU(ctx, position, box, radius, props)
    }
}

export async function GaussianDensityCPU(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps): Promise<DensityData> {
    const { resolution, radiusOffset, smoothness } = props

    const { indices, x, y, z } = position
    const n = OrderedSet.size(indices)

    const v = Vec3.zero()
    const p = Vec3.zero()

    const pad = (radiusOffset + 3) * 3 // TODO calculate max radius
    const expandedBox = Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad))
    const extent = Vec3.sub(Vec3.zero(), expandedBox.max, expandedBox.min)
    const min = expandedBox.min

    const delta = getDelta(Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad)), resolution)
    const dim = Vec3.zero()
    Vec3.ceil(dim, Vec3.mul(dim, extent, delta))

    const space = Tensor.Space(dim, [0, 1, 2], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    const idData = space.create()
    const idField = Tensor.create(space, idData)

    const densData = space.create()

    const c = Vec3.zero()

    const alpha = smoothness

    const _r2 = (radiusOffset + 1.4 * 2)
    const _radius2 = Vec3.create(_r2, _r2, _r2)
    Vec3.mul(_radius2, _radius2, delta)
    const updateChunk = Math.ceil(10000 / (_radius2[0] * _radius2[1] * _radius2[2]))

    const beg = Vec3.zero()
    const end = Vec3.zero()

    const gridPad = 1 / Math.max(...delta)

    console.time('gaussian density cpu')
    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i)

        Vec3.set(v, x[j], y[j], z[j])

        Vec3.sub(v, v, min)
        Vec3.mul(c, v, delta)

        const rad = radius(j) + radiusOffset
        const rSq = rad * rad

        const r2 = radiusOffset + rad * 2 + gridPad
        const rad2 = Vec3.create(r2, r2, r2)
        Vec3.mul(rad2, rad2, delta)
        const r2sq = r2 * r2

        const [ begX, begY, begZ ] = Vec3.floor(beg, Vec3.sub(beg, c, rad2))
        const [ endX, endY, endZ ] = Vec3.ceil(end, Vec3.add(end, c, rad2))

        for (let xi = begX; xi < endX; ++xi) {
            for (let yi = begY; yi < endY; ++yi) {
                for (let zi = begZ; zi < endZ; ++zi) {
                    Vec3.set(p, xi, yi, zi)
                    Vec3.div(p, p, delta)
                    const distSq = Vec3.squaredDistance(p, v)
                    if (distSq <= r2sq) {
                        const dens = Math.exp(-alpha * (distSq / rSq))
                        space.add(data, xi, yi, zi, dens)
                        if (dens > space.get(densData, xi, yi, zi)) {
                            space.set(densData, xi, yi, zi, dens)
                            space.set(idData, xi, yi, zi, i)
                        }
                    }
                }
            }
        }

        if (i % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'filling density grid', current: i, max: n })
        }
    }
    console.timeEnd('gaussian density cpu')

    const transform = Mat4.identity()
    Mat4.fromScaling(transform, Vec3.inverse(Vec3.zero(), delta))
    Mat4.setTranslation(transform, expandedBox.min)

    return { field, idField, transform }
}

export async function GaussianDensityGPU(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps) { // }: Promise<DensityData> {
    const { resolution, radiusOffset, smoothness, readSlices } = props

    const { indices, x, y, z } = position
    const n = OrderedSet.size(indices)

    const positions = new Float32Array(n * 3)
    const radii = new Float32Array(n)

    const pad = (radiusOffset + 3) * 3 // TODO calculate max radius
    const expandedBox = Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad));
    const extent = Vec3.sub(Vec3.zero(), expandedBox.max, expandedBox.min)

    const delta = getDelta(Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad)), resolution)
    const dim = Vec3.zero()
    Vec3.ceil(dim, Vec3.mul(dim, extent, delta))

    const _r2 = (radiusOffset + 1.4 * 2)
    const _radius2 = Vec3.create(_r2, _r2, _r2)
    Vec3.mul(_radius2, _radius2, delta)
    const updateChunk = Math.ceil(10000 / (_radius2[0] * _radius2[1] * _radius2[2]))

    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i);

        positions[i * 3] = x[j]
        positions[i * 3 + 1] = y[j]
        positions[i * 3 + 2] = z[j]
        radii[i] = radius(j) + radiusOffset

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'preparing density data', current: i, max: n })
        }
    }

    //

    const values: GaussianDensityValues = {
        drawCount: ValueCell.create(n),
        instanceCount: ValueCell.create(1),

        aRadius: ValueCell.create(radii),
        aPosition: ValueCell.create(positions),

        uCurrentSlice: ValueCell.create(0),
        uCurrentX: ValueCell.create(0),
        uCurrentY: ValueCell.create(0),
        uBboxMin: ValueCell.create(expandedBox.min),
        uBboxMax: ValueCell.create(expandedBox.max),
        uBboxSize: ValueCell.create(extent),
        uGridDim: ValueCell.create(dim),
        uAlpha: ValueCell.create(smoothness),
    }
    const state: RenderableState = {
        visible: true,
        depthMask: false
    }

    // TODO do in OffscreenCanvas (https://www.chromestatus.com/feature/5681560598609920)
    const webgl = getWebGLContext()

    const renderObject = createGaussianDensityRenderObject(values, state)
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

    //

    const space = Tensor.Space(dim, [2, 1, 0], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    const idData = space.create()
    const idField = Tensor.create(space, idData)

    //

    const { gl } = webgl

    const program = renderable.getProgram('draw')
    const renderTarget = createRenderTarget(webgl, fboTexDimX, fboTexDimY)

    program.use()
    renderTarget.bind()

    gl.disable(gl.CULL_FACE)
    gl.frontFace(gl.CCW)
    gl.cullFace(gl.BACK)

    gl.depthMask(true)
    gl.clearColor(0, 0, 0, 1)
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.depthMask(false)

    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
    gl.blendEquation(gl.FUNC_ADD)
    gl.enable(gl.BLEND)

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
        ValueCell.update(values.uCurrentSlice, i)
        ValueCell.update(values.uCurrentX, currX)
        ValueCell.update(values.uCurrentY, currY)
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

        if (i % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'filling density grid', current: i, max: n })
        }
    }
    console.timeEnd('gpu gaussian density slices')

    //

    if (!readSlices) {
        console.time('gpu gaussian density full')
        renderTarget.getBuffer()
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
                    data[idx] = renderTarget.image.array[4 * (tmpCol * dim[0] + (iy + tmpRow) * fboTexDimX + ix)] / 255
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

    return { field, idField, transform }
}

let webglContext: Context
function getWebGLContext() {
    if (webglContext) return webglContext
    const canvas = document.createElement('canvas')
    const gl = canvas.getContext('webgl', {
        alpha: false,
        antialias: true,
        depth: true,
        preserveDrawingBuffer: true
    })
    if (!gl) throw new Error('Could not create a WebGL rendering context')
    webglContext = createContext(gl)
    return webglContext
}