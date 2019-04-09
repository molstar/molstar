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
    const r = state.position.radius[j] + state.probeRadius
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
    const { position, probeRadius, lookup3d, min, delta, space, data, idData } = state

    const { indices, x, y, z, radius } = position
    const n = OrderedSet.size(indices)

    const v = Vec3()
    const p = Vec3()
    const c = Vec3()
    const sp = Vec3()

    const beg = Vec3()
    const end = Vec3()

    // const gridPad = 1 / Math.max(...delta)

    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i)

        Vec3.set(v, x[j], y[j], z[j])

        state.neighbours = lookup3d.find(v[0], v[1], v[2], radius[j] + probeRadius)

        Vec3.sub(v, v, min)
        Vec3.mul(c, v, delta)

        const rad = radius[j] + probeRadius
        const rSq = rad * rad

        const r2 = rad // * 2 + gridPad
        const rad2 = Vec3.create(r2, r2, r2)
        Vec3.mul(rad2, rad2, delta)

        const [ begX, begY, begZ ] = Vec3.floor(beg, Vec3.sub(beg, c, rad2))
        const [ endX, endY, endZ ] = Vec3.ceil(end, Vec3.add(end, c, rad2))

        for (let xi = begX; xi < endX; ++xi) {
            for (let yi = begY; yi < endY; ++yi) {
                for (let zi = begZ; zi < endZ; ++zi) {
                    Vec3.set(p, xi, yi, zi)
                    Vec3.div(p, p, delta)

                    const dv = Vec3()
                    Vec3.sub(dv, p, v)

                    // const distSq = Vec3.squaredDistance(p, v)
                    const dSq = Vec3.squaredMagnitude(dv)

                    if (dSq < rSq) {
                        const val = space.get(data, xi, yi, zi)
                        if (val < 0.0) {
                            // Unvisited, make positive
                            space.set(data, xi, yi, zi, -val)
                        }

                        // Project on to the surface of the sphere
                        // sp is the projected point ( dx, dy, dz ) * ( ra / d )
                        // const dist = Math.sqrt(distSq)
                        const d = Math.sqrt(dSq)
                        const ap = rad / d
                        Vec3.scale(sp, dv, ap)
                        Vec3.add(sp, sp, v)
                        Vec3.add(sp, sp, min)
                        // Vec3.add(sp, v, Vec3.setMagnitude(sp, Vec3.sub(sp, p, v), rad - probeRadius))

                        if (obscured(state, sp[0], sp[1], sp[2], i, -1) === -1) {
                            const dd = rad - d
                            if (dd < space.get(data, xi, yi, zi)) {
                                space.set(data, xi, yi, zi, dd)
                                space.set(idData, xi, yi, zi, i)
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
const v = Vec3()
const p = Vec3()
const c = Vec3()
const beg = Vec3()
const end = Vec3()
const radVec = Vec3()
function projectTorus (state: MolSurfCalcState, a: number, b: number) {
    const { position, min, delta, space, data, idData } = state
    const { cosTable, sinTable, probePositions, probeRadius, resolution } = state

    const r1 = position.radius[a] + probeRadius
    const r2 = position.radius[b] + probeRadius
    const dx = atob[0] = position.x[b] - position.x[a]
    const dy = atob[1] = position.y[b] - position.y[a]
    const dz = atob[2] = position.z[b] - position.z[a]
    const dSq = dx * dx + dy * dy + dz * dz

    // This check now redundant as already done in AVHash.withinRadii
    if (dSq > ((r1 + r2) * (r1 + r2))) { return }

    const d = Math.sqrt(dSq)

    // Find angle between a->b vector and the circle
    // of their intersection by cosine rule
    const cosA = (r1 * r1 + d * d - r2 * r2) / (2.0 * r1 * d)

    // distance along a->b at intersection
    const dmp = r1 * cosA

    Vec3.normalize(atob, atob)

    // Create normal to line
    normalToLine(n1, atob)
    Vec3.normalize(n1, n1)

    // Cross together for second normal vector
    Vec3.cross(n2, atob, n1)
    Vec3.normalize(n2, n2)

    // r is radius of circle of intersection
    const rInt = Math.sqrt(r1 * r1 - dmp * dmp)

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

            Vec3.set(v, px, py, pz)

            Vec3.sub(v, v, min)
            Vec3.mul(c, v, delta)

            const rad = probeRadius / resolution
            Vec3.set(radVec, rad, rad, rad)
            Vec3.mul(radVec, radVec, delta)

            const [ begX, begY, begZ ] = Vec3.floor(beg, Vec3.sub(beg, c, radVec))
            const [ endX, endY, endZ ] = Vec3.ceil(end, Vec3.add(end, c, radVec))

            for (let xi = begX; xi < endX; ++xi) {
                for (let yi = begY; yi < endY; ++yi) {
                    for (let zi = begZ; zi < endZ; ++zi) {
                        Vec3.set(p, xi, yi, zi)
                        Vec3.div(p, p, delta)

                        const dv = Vec3()
                        Vec3.sub(dv, v, p)
                        const dSq = Vec3.squaredMagnitude(dv)

                        // const distSq = Vec3.squaredDistance(p, v)
                        const current = space.get(data, xi, yi, zi)

                        if (current > 0.0 && dSq < (current * current)) {
                            space.set(data, xi, yi, zi, Math.sqrt(dSq))
                            // Is this grid point closer to a or b?
                            // Take dot product of atob and gridpoint->p (dx, dy, dz)
                            const dp = dx * atob[0] + dy * atob [1] + dz * atob[2]
                            space.set(idData, xi, yi, zi, dp < 0.0 ? b : a)
                        }
                    }
                }
            }
        }
    }
}

async function projectTorii (ctx: RuntimeContext, state: MolSurfCalcState) {
    const { n, lookup3d, position, probeRadius, resolution } = state
    const { x, y, z, radius } = position
    for (let i = 0; i < n; ++i) {
        state.neighbours = lookup3d.find(x[i], y[i], z[i], radius[i] + probeRadius / resolution)
        for (let j = 0, jl = state.neighbours.count; j < jl; ++j) {
            const ib = state.neighbours.indices[j]
            if (i < ib) projectTorus(state, i, ib)
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
    min: Vec3

    maxRadius: number

    n: number
    resolution: number
    probeRadius: number

    /** Angle lookup tables */
    cosTable: Float32Array
    sinTable: Float32Array
    probePositions: number

    expandedBox: Box3D
    space: Tensor.Space
    data: Tensor.Data
    idData: Tensor.Data
}

async function createState(ctx: RuntimeContext, position: Required<PositionData>, maxRadius: number, props: MolecularSurfaceCalculationProps): Promise<MolSurfCalcState> {
    const { resolution, probeRadius, probePositions } = props

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

    const { cosTable, sinTable } = getAngleTables(probePositions)

    const space = Tensor.Space(dim, [0, 1, 2], Float32Array)
    const data = space.create()
    const idData = space.create()

    fillUniform(data, -1001.0)
    fillUniform(idData, -1)

    return {
        lastClip: -1,
        neighbours: lookup3d.find(0, 0, 0, 0),

        lookup3d,
        position,
        delta,
        min,

        maxRadius,

        n,
        resolution,
        probeRadius,

        cosTable,
        sinTable,
        probePositions,

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