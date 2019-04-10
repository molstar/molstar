/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D, fillGridDim } from '../../geometry';
import { Vec3, Mat4, Tensor } from '../../linear-algebra';
import { RuntimeContext } from 'mol-task';
import { PositionData, DensityData } from '../common';
import { OrderedSet } from 'mol-data/int';
import { GaussianDensityProps, getDelta } from '../gaussian-density';

export async function GaussianDensityCPU(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps): Promise<DensityData> {
    const { resolution, radiusOffset, smoothness } = props
    const scaleFactor = 1 / resolution

    const { indices, x, y, z } = position
    const n = OrderedSet.size(indices)
    const radii = new Float32Array(n)

    let maxRadius = 0
    for (let i = 0; i < n; ++i) {
        const r = radius(OrderedSet.getAt(indices, i)) + radiusOffset
        if (maxRadius < r) maxRadius = r
        radii[i] = r

        if (i % 100000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'calculating max radius', current: i, max: n })
        }
    }

    const pad = maxRadius * 2 + resolution
    const expandedBox = Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad))
    const extent = Vec3.sub(Vec3(), expandedBox.max, expandedBox.min)
    const min = expandedBox.min

    const delta = getDelta(Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad)), resolution)
    const dim = Vec3()
    Vec3.ceil(dim, Vec3.mul(dim, extent, delta))
    // console.log('grid dim cpu', dim)

    const space = Tensor.Space(dim, [0, 1, 2], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    const idData = space.create()
    const idField = Tensor.create(space, idData)

    const [ dimX, dimY, dimZ ] = dim
    const iu = dimZ, iv = dimY, iuv = iu * iv

    const gridx = fillGridDim(dim[0], min[0], resolution)
    const gridy = fillGridDim(dim[1], min[1], resolution)
    const gridz = fillGridDim(dim[2], min[2], resolution)

    const densData = space.create()

    const alpha = smoothness
    const updateChunk = Math.ceil(1000000 / (Math.pow(maxRadius * 2, 3) * resolution))

    // console.time('gaussian density cpu')
    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i)
        const vx = x[j], vy = y[j], vz = z[j]

        const rad = radii[i]
        const rSq = rad * rad
        const rSqInv = 1 / rSq

        const r2 = rad * 2
        const r2sq = r2 * r2

        // Number of grid points, round this up...
        const ng = Math.ceil(r2 * scaleFactor)

        // Center of the atom, mapped to grid points (take floor)
        const iax = Math.floor(scaleFactor * (vx - min[0]))
        const iay = Math.floor(scaleFactor * (vy - min[1]))
        const iaz = Math.floor(scaleFactor * (vz - min[2]))

        // Extents of grid to consider for this atom
        const begX = Math.max(0, iax - ng)
        const begY = Math.max(0, iay - ng)
        const begZ = Math.max(0, iaz - ng)

        // Add two to these points:
        // - iax are floor'd values so this ensures coverage
        // - these are loop limits (exclusive)
        const endX = Math.min(dimX, iax + ng + 2)
        const endY = Math.min(dimY, iay + ng + 2)
        const endZ = Math.min(dimZ, iaz + ng + 2)

        for (let xi = begX; xi < endX; ++xi) {
            const dx = gridx[xi] - vx
            const xIdx = xi * iuv
            for (let yi = begY; yi < endY; ++yi) {
                const dy = gridy[yi] - vy
                const dxySq = dx * dx + dy * dy
                const xyIdx = yi * iu + xIdx
                for (let zi = begZ; zi < endZ; ++zi) {
                    const dz = gridz[zi] - vz
                    const dSq = dxySq + dz * dz
                    if (dSq <= r2sq) {
                        const dens = Math.exp(-alpha * (dSq * rSqInv))
                        const idx = zi + xyIdx
                        data[idx] += dens
                        if (dens > densData[idx]) {
                            densData[idx] = dens
                            idData[idx] = i
                        }
                    }
                }
            }
        }

        if (i % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'filling density grid', current: i, max: n })
        }
    }
    // console.timeEnd('gaussian density cpu')

    const transform = Mat4.identity()
    Mat4.fromScaling(transform, Vec3.inverse(Vec3.zero(), delta))
    Mat4.setTranslation(transform, expandedBox.min)

    return { field, idField, transform }
}