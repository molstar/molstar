/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Fred Ludlow <fred.ludlow@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * ported from NGL (https://github.com/arose/ngl), licensed under MIT
 */

import { fillUniform } from 'mol-util/array';
import { Vec3, Tensor } from 'mol-math/linear-algebra';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Lookup3D, Result } from './lookup3d/common';
import { RuntimeContext } from 'mol-task';
import { OrderedSet } from 'mol-data/int';
import { PositionData } from './common';
import { Mat4 } from 'mol-math/linear-algebra/3d';
import { Box3D, GridLookup3D } from 'mol-math/geometry';
import { getDelta } from './gaussian-density';

function normalToLine (out: Vec3, p: Vec3) {
    out[0] = out[1] = out[2] = 1.0
    if (p[0] !== 0) {
        out[0] = (p[1] + p[2]) / -p[0]
    } else if (p[1] !== 0) {
        out[1] = (p[0] + p[2]) / -p[1]
    } else if (p[2] !== 0) {
        out[2] = (p[0] + p[1]) / -p[2]
    }
    return out
}

type AnglesTables = { cosTable: Float32Array, sinTable: Float32Array }
function getAngleTables (probePositions: number): AnglesTables {
    let theta = 0.0
    const step = 2 * Math.PI / probePositions

    const cosTable = new Float32Array(probePositions)
    const sinTable = new Float32Array(probePositions)
    for (let i = 0; i < probePositions; i++) {
        cosTable[i] = Math.cos(theta)
        sinTable[i] = Math.sin(theta)
        theta += step
    }
    return { cosTable, sinTable}
}

function fillGridDim(a: Float32Array, start: number, step: number) {
    for (let i = 0; i < a.length; i++) {
        a[i] = start + (step * i)
    }
}

/**
 * Is the point at x,y,z obscured by any of the atoms specifeid by indices in neighbours.
 * Ignore indices a and b (these are the relevant atoms in projectPoints/Torii)
 *
 * Cache the last clipped atom (as very often the same one in subsequent calls)
 */
function obscured (state: MolSurfCalcState, x: number, y: number, z: number, a: number, b: number) {
    if (state.lastClip !== -1) {
        const ai = state.lastClip
        if (ai !== a && ai !== b && singleAtomObscures(state, ai, x, y, z)) {
            return ai
        } else {
            state.lastClip = -1
        }
    }

    for (let j = 0, jl = state.neighbours.count; j < jl; ++j) {
        const ai = state.neighbours.indices[j]
        if (ai !== a && ai !== b && singleAtomObscures(state, ai, x, y, z)) {
            state.lastClip = ai
            return ai
        }
    }

    return -1
}

function singleAtomObscures (state: MolSurfCalcState, ai: number, x: number, y: number, z: number) {
    const j = OrderedSet.getAt(state.position.indices, ai)
    const r = state.position.radius[j]
    const dx = state.position.x[j] - x
    const dy = state.position.y[j] - y
    const dz = state.position.z[j] - z
    const dSq = dx * dx + dy * dy + dz * dz
    return dSq < (r * r)
}

/**
 * For each atom:
 *     Iterate over a subsection of the grid, for each point:
 *         If current value < 0.0, unvisited, set positive
 *
 *         In any case: Project this point onto surface of the atomic sphere
 *         If this projected point is not obscured by any other atom
 *             Calculate delta distance and set grid value to minimum of
 *             itself and delta
 */
async function projectPoints (ctx:  RuntimeContext, state: MolSurfCalcState) {
    const { position, lookup3d, min, space, data, idData, scaleFactor } = state
    const { gridx, gridy, gridz } = state

    const { indices, x, y, z, radius } = position
    const n = OrderedSet.size(indices)

    const [ dimX, dimY, dimZ ] = space.dimensions
    const iu = dimZ, iv = dimY, iuv = iu * iv

    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i)
        const vx = x[j], vy = y[j], vz = z[j]
        const rad = radius[j]
        const rSq = rad * rad

        state.neighbours = lookup3d.find(vx, vy, vz, rad)

        // Number of grid points, round this up...
        const ng = Math.ceil(rad * scaleFactor)

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

                    if (dSq < rSq) {
                        const idx = zi + xyIdx

                        // if unvisited, make positive
                        if (data[idx] < 0.0) data[idx] *= -1

                        // Project on to the surface of the sphere
                        // sp is the projected point ( dx, dy, dz ) * ( ra / d )
                        const d = Math.sqrt(dSq)
                        const ap = rad / d
                        const spx = dx * ap + vx
                        const spy = dy * ap + vy
                        const spz = dz * ap + vz

                        if (obscured(state, spx, spy, spz, i, -1) === -1) {
                            const dd = rad - d
                            if (dd < data[idx]) {
                                data[idx] = dd
                                idData[idx] = i
                            }
                        }
                    }
                }
            }
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'projecting points', current: i, max: n })
        }
    }
}

// Vectors for Torus Projection
const atob = Vec3()
const mid = Vec3()
const n1 = Vec3()
const n2 = Vec3()
function projectTorus (state: MolSurfCalcState, a: number, b: number) {
    const { position, min, space, data, idData } = state
    const { cosTable, sinTable, probePositions, probeRadius, scaleFactor } = state
    const { gridx, gridy, gridz } = state

    const [ dimX, dimY, dimZ ] = space.dimensions
    const iu = dimZ, iv = dimY, iuv = iu * iv

    const ng = Math.max(5, 2 + Math.floor(probeRadius * scaleFactor))

    const rA = position.radius[a]
    const rB = position.radius[b]
    const dx = atob[0] = position.x[b] - position.x[a]
    const dy = atob[1] = position.y[b] - position.y[a]
    const dz = atob[2] = position.z[b] - position.z[a]
    const dSq = dx * dx + dy * dy + dz * dz

    // This check now redundant as already done in AVHash.withinRadii
    if (dSq > ((rA + rB) * (rA + rB))) { return }

    const d = Math.sqrt(dSq)

    // Find angle between a->b vector and the circle
    // of their intersection by cosine rule
    const cosA = (rA * rA + d * d - rB * rB) / (2.0 * rA * d)

    // distance along a->b at intersection
    const dmp = rA * cosA

    Vec3.normalize(atob, atob)

    // Create normal to line
    normalToLine(n1, atob)
    Vec3.normalize(n1, n1)

    // Cross together for second normal vector
    Vec3.cross(n2, atob, n1)
    Vec3.normalize(n2, n2)

    // r is radius of circle of intersection
    const rInt = Math.sqrt(rA * rA - dmp * dmp)

    Vec3.scale(n1, n1, rInt)
    Vec3.scale(n2, n2, rInt)
    Vec3.scale(atob, atob, dmp)

    mid[0] = atob[0] + position.x[a]
    mid[1] = atob[1] + position.y[a]
    mid[2] = atob[2] + position.z[a]

    state.lastClip = -1

    for (let i = 0; i < probePositions; ++i) {
        const cost = cosTable[i]
        const sint = sinTable[i]

        const px = mid[0] + cost * n1[0] + sint * n2[0]
        const py = mid[1] + cost * n1[1] + sint * n2[1]
        const pz = mid[2] + cost * n1[2] + sint * n2[2]

        if (obscured(state, px, py, pz, a, b) === -1) {
            const iax = Math.floor(scaleFactor * (px - min[0]))
            const iay = Math.floor(scaleFactor * (py - min[1]))
            const iaz = Math.floor(scaleFactor * (pz - min[2]))

            const begX = Math.max(0, iax - ng)
            const begY = Math.max(0, iay - ng)
            const begZ = Math.max(0, iaz - ng)

            const endX = Math.min(dimX, iax + ng + 2)
            const endY = Math.min(dimY, iay + ng + 2)
            const endZ = Math.min(dimZ, iaz + ng + 2)

            for (let xi = begX; xi < endX; ++xi) {
                const dx = px - gridx[xi]
                const xIdx = xi * iuv

                for (let yi = begY; yi < endY; ++yi) {
                    const dy = py - gridy[yi]
                    const dxySq = dx * dx + dy * dy
                    const xyIdx = yi * iu + xIdx

                    for (let zi = begZ; zi < endZ; ++zi) {
                        const dz = pz - gridz[zi]
                        const dSq = dxySq + dz * dz

                        const idx = zi + xyIdx
                        const current = data[idx]

                        if (current > 0.0 && dSq < (current * current)) {
                            data[idx] = Math.sqrt(dSq)
                            // Is this grid point closer to a or b?
                            // Take dot product of atob and gridpoint->p (dx, dy, dz)
                            const dp = dx * atob[0] + dy * atob[1] + dz * atob[2]
                            idData[idx] = dp < 0.0 ? b : a
                        }
                    }
                }
            }
        }
    }
}

async function projectTorii (ctx: RuntimeContext, state: MolSurfCalcState) {
    const { n, lookup3d, position } = state
    const { x, y, z, indices, radius } = position
    for (let i = 0; i < n; ++i) {
        const k = OrderedSet.getAt(indices, i)
        state.neighbours = lookup3d.find(x[k], y[k], z[k], radius[k])
        for (let j = 0, jl = state.neighbours.count; j < jl; ++j) {
            const l = state.neighbours.indices[j]
            if (k < l) projectTorus(state, k, l)
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'projecting torii', current: i, max: n })
        }
    }
}

//

interface MolSurfCalcState {
    /** Cached last value for obscured test */
    lastClip: number
    /** Neighbours as transient result array from lookup3d */
    neighbours: Result<number>

    lookup3d: Lookup3D
    position: Required<PositionData>
    delta: Vec3
    invDelta: Vec3
    min: Vec3

    maxRadius: number

    n: number
    resolution: number
    scaleFactor: number
    probeRadius: number

    /** Angle lookup tables */
    cosTable: Float32Array
    sinTable: Float32Array
    probePositions: number

    /** grid lookup tables */
    gridx: Float32Array
    gridy: Float32Array
    gridz: Float32Array

    expandedBox: Box3D
    space: Tensor.Space
    data: Tensor.Data
    idData: Tensor.Data
}

async function createState(ctx: RuntimeContext, position: Required<PositionData>, maxRadius: number, props: MolecularSurfaceCalculationProps): Promise<MolSurfCalcState> {
    const { resolution, probeRadius, probePositions } = props

    const scaleFactor = 1 / resolution

    const lookup3d = GridLookup3D(position)
    const box = lookup3d.boundary.box
    const { indices } = position
    const n = OrderedSet.size(indices)

    const pad = maxRadius * 2 + resolution
    const expandedBox = Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad))
    const extent = Vec3.sub(Vec3.zero(), expandedBox.max, expandedBox.min)
    const min = expandedBox.min

    const delta = getDelta(Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad)), resolution)
    const dim = Vec3.zero()
    Vec3.ceil(dim, Vec3.mul(dim, extent, delta))
    console.log('grid dim surf', dim)
    const invDelta = Vec3.inverse(Vec3(), delta)

    const { cosTable, sinTable } = getAngleTables(probePositions)

    const space = Tensor.Space(dim, [0, 1, 2], Float32Array)
    const data = space.create()
    const idData = space.create()

    fillUniform(data, -1001.0)
    fillUniform(idData, -1)

    const gridx = new Float32Array(dim[0])
    const gridy = new Float32Array(dim[1])
    const gridz = new Float32Array(dim[2])

    fillGridDim(gridx, min[0], resolution)
    fillGridDim(gridy, min[1], resolution)
    fillGridDim(gridz, min[2], resolution)

    return {
        lastClip: -1,
        neighbours: lookup3d.find(0, 0, 0, 0),

        lookup3d,
        position,
        delta,
        invDelta,
        min,

        maxRadius,

        n,
        resolution,
        scaleFactor,
        probeRadius,

        cosTable,
        sinTable,
        probePositions,

        gridx,
        gridy,
        gridz,

        expandedBox,
        space,
        data,
        idData,
    }
}

//

export const MolecularSurfaceCalculationParams = {
    resolution: PD.Numeric(0.5, { min: 0.01, max: 10, step: 0.01 }),
    probeRadius: PD.Numeric(1.4, { min: 0, max: 10, step: 0.1 }),
    probePositions: PD.Numeric(30, { min: 12, max: 90, step: 1 }),
}
export const DefaultMolecularSurfaceCalculationProps = PD.getDefaultValues(MolecularSurfaceCalculationParams)
export type MolecularSurfaceCalculationProps = typeof DefaultMolecularSurfaceCalculationProps


export async function calcMolecularSurface(ctx: RuntimeContext, position: Required<PositionData>, maxRadius: number,  props: MolecularSurfaceCalculationProps) {
    // Field generation method adapted from AstexViewer (Mike Hartshorn) by Fred Ludlow.
    // Other parts based heavily on NGL (Alexander Rose) EDT Surface class

    console.time('MolecularSurface')

    console.time('MolecularSurface createState')
    const state = await createState(ctx, position, maxRadius, props)
    console.timeEnd('MolecularSurface createState')

    console.time('MolecularSurface projectPoints')
    await projectPoints(ctx, state)
    console.timeEnd('MolecularSurface projectPoints')

    console.time('MolecularSurface projectTorii')
    await projectTorii(ctx, state)
    console.timeEnd('MolecularSurface projectTorii')

    console.timeEnd('MolecularSurface')

    const field = Tensor.create(state.space, state.data)
    const idField = Tensor.create(state.space, state.idData)

    const transform = Mat4.identity()
    Mat4.fromScaling(transform, Vec3.inverse(Vec3.zero(), state.delta))
    Mat4.setTranslation(transform, state.expandedBox.min)
    console.log({ field, idField, transform, state })
    return { field, idField, transform }
}