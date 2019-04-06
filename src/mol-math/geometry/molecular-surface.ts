/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Fred Ludlow <fred.ludlow@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * ported from NGL (https://github.com/arose/ngl), licensed under MIT
 */

import { fillUniform } from 'mol-util/array';
import { Vec3 } from 'mol-math/linear-algebra';
import { NumberArray } from 'mol-util/type-helpers';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { number } from 'prop-types';
import { Lookup3D, Result } from './lookup3d/common';
import { RuntimeContext } from 'mol-task';
import { OrderedSet } from 'mol-data/int';
import { PositionData } from './common';

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

function fillGridDim (a: Float32Array, start: number, step: number) {
    for (let i = 0; i < a.length; i++) {
        a[i] = start + (step * i)
    }
}

type AnglesTables = { cosTable: Float32Array, sinTable: Float32Array }
function getAngleTables (probePositions: number): AnglesTables {
    let theta = 0.0
    const step = 2 * Math.PI / probePositions

    const cosTable = new Float32Array(probePositions)
    const sinTable = new Float32Array(probePositions)
    for (let i = 0; i < probePositions; i++) {
        cosTable[ i ] = Math.cos(theta)
        sinTable[ i ] = Math.sin(theta)
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
    let ai: number

    if (state.lastClip !== -1) {
        ai = state.lastClip
        if (ai !== a && ai !== b && singleAtomObscures(state, ai, x, y, z)) {
            return ai
        } else {
            state.lastClip = -1
        }
    }

    let ni = 0
    ai = state.neighbours[ni]
    while (ai >= 0) {
        if (ai !== a && ai !== b && singleAtomObscures(state, ai, x, y, z)) {
            state.lastClip = ai
            return ai
        }
        ai = state.neighbours[++ni]
    }

    state.lastClip = -1

    return -1
}

function singleAtomObscures (state: MolSurfCalcState, ai: number, x: number, y: number, z: number) {
    const ra2 = state.radiusSq[ai]
    const dx = state.xCoord[ai] - x
    const dy = state.yCoord[ai] - y
    const dz = state.zCoord[ai] - z
    const d2 = dx * dx + dy * dy + dz * dz
    return d2 < ra2
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
    const { position, radius, scaleFactor, lookup3D, min, delta } = state

    const { indices, x, y, z } = position
    const n = OrderedSet.size(indices)

    const v = Vec3()
    const p = Vec3()
    const c = Vec3()

    const beg = Vec3()
    const end = Vec3()

    const gridPad = 1 / Math.max(...delta)

    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i)

        Vec3.set(v, x[j], y[j], z[j])

        Vec3.sub(v, v, min)
        Vec3.mul(c, v, delta)

        const rad = radius(j)
        const rSq = rad * rad

        const r2 = rad * 2 + gridPad
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

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'projecting points', current: i, max: n })
        }
    }

    for (let i = 0; i < nAtoms; i++) {
        const ax = xCoord[ i ]
        const ay = yCoord[ i ]
        const az = zCoord[ i ]
        const ar = radius[ i ]
        const ar2 = radiusSq[ i ]

        state.neighbours = lookup3D.find(ax, ay, az, ar)

        // Number of grid points, round this up...
        const ng = Math.ceil(ar * scaleFactor)

        // Center of the atom, mapped to grid points (take floor)
        const iax = Math.floor(scaleFactor * (ax - min[ 0 ]))
        const iay = Math.floor(scaleFactor * (ay - min[ 1 ]))
        const iaz = Math.floor(scaleFactor * (az - min[ 2 ]))

        // Extents of grid to consider for this atom
        const minx = Math.max(0, iax - ng)
        const miny = Math.max(0, iay - ng)
        const minz = Math.max(0, iaz - ng)

        // Add two to these points:
        // - iax are floor'd values so this ensures coverage
        // - these are loop limits (exclusive)
        const maxx = Math.min(dim[ 0 ], iax + ng + 2)
        const maxy = Math.min(dim[ 1 ], iay + ng + 2)
        const maxz = Math.min(dim[ 2 ], iaz + ng + 2)

        for (let ix = minx; ix < maxx; ix++) {
            const dx = gridx[ ix ] - ax
            const xoffset = dim[ 1 ] * dim[ 2 ] * ix

            for (let iy = miny; iy < maxy; iy++) {
                const dy = gridy[ iy ] - ay
                const dxy2 = dx * dx + dy * dy
                const xyoffset = xoffset + dim[ 2 ] * iy

                for (let iz = minz; iz < maxz; iz++) {
                    const dz = gridz[ iz ] - az
                    const d2 = dxy2 + dz * dz

                    if (d2 < ar2) {
                        const idx = iz + xyoffset

                        if (grid[idx] < 0.0) {
                            // Unvisited, make positive
                            grid[ idx ] = -grid[ idx ]
                        }
                        // Project on to the surface of the sphere
                        // sp is the projected point ( dx, dy, dz ) * ( ra / d )
                        const d = Math.sqrt(d2)
                        const ap = ar / d
                        let spx = dx * ap
                        let spy = dy * ap
                        let spz = dz * ap

                        spx += ax
                        spy += ay
                        spz += az

                        if (obscured(state, spx, spy, spz, i, -1) === -1) {
                            const dd = ar - d
                            if (dd < grid[ idx ]) {
                                grid[ idx ] = dd
                                atomIndex[ idx ] = i
                            }
                        }
                    }
                }
            }
        }
    }
}

// Vectors for Torus Projection
const atob = Vec3()
const mid = Vec3()
const n1 = Vec3()
const n2 = Vec3()
function projectTorus (state: MolSurfCalcState, a: number, b: number) {
    const r1 = state.radius[a]
    const r2 = state.radius[b]
    const dx = atob[0] = state.xCoord[b] - state.xCoord[a]
    const dy = atob[1] = state.yCoord[b] - state.yCoord[a]
    const dz = atob[2] = state.zCoord[b] - state.zCoord[a]
    const d2 = dx * dx + dy * dy + dz * dz

    // This check now redundant as already done in AVHash.withinRadii
    // if (d2 > ((r1 + r2) * (r1 + r2))){ return; }

    const d = Math.sqrt(d2)

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

    mid[0] = atob[0] + state.xCoord[a]
    mid[1] = atob[1] + state.yCoord[a]
    mid[2] = atob[2] + state.zCoord[a]

    state.lastClip = -1

    const { ngTorus, cosTable, sinTable, scaleFactor } = state

    for (let i = 0; i < state.probePositions; i++) {
        const cost = cosTable[ i ]
        const sint = sinTable[ i ]

        const px = mid[0] + cost * n1[0] + sint * n2[0]
        const py = mid[1] + cost * n1[1] + sint * n2[1]
        const pz = mid[2] + cost * n1[2] + sint * n2[2]

        if (obscured(state, px, py, pz, a, b) === -1) {
            // As above, iterate over our grid...
            // px, py, pz in grid coords
            const iax = Math.floor(scaleFactor * (px - min[0]))
            const iay = Math.floor(scaleFactor * (py - min[1]))
            const iaz = Math.floor(scaleFactor * (pz - min[2]))

            const minx = Math.max(0, iax - ngTorus)
            const miny = Math.max(0, iay - ngTorus)
            const minz = Math.max(0, iaz - ngTorus)

            const maxx = Math.min(dim[0], iax + ngTorus + 2)
            const maxy = Math.min(dim[1], iay + ngTorus + 2)
            const maxz = Math.min(dim[2], iaz + ngTorus + 2)

            for (let ix = minx; ix < maxx; ix++) {
                const dx = px - gridx[ ix ]
                const xoffset = dim[1] * dim[2] * ix

                for (let iy = miny; iy < maxy; iy++) {
                    const dy = py - gridy[iy]
                    const  dxy2 = dx * dx + dy * dy
                    const  xyoffset = xoffset + dim[2] * iy

                    for (let iz = minz; iz < maxz; iz++) {
                        const dz = pz - gridz[iz]
                        const d2 = dxy2 + dz * dz
                        const  idx = iz + xyoffset
                        const  current = grid[idx]

                        if (current > 0.0 && d2 < (current * current)) {
                            grid[idx] = Math.sqrt(d2)
                            // Is this grid point closer to a or b?
                            // Take dot product of atob and gridpoint->p (dx, dy, dz)
                            const dp = dx * atob[0] + dy * atob [1] + dz * atob[2]
                            atomIndex[idx] = dp < 0.0 ? b : a
                        }
                    }
                }
            }
        }
    }
}

function projectTorii (state: MolSurfCalcState) {
    const { n: nAtoms, neighbours, hash, xCoord, yCoord, zCoord, radius } = state
    for (let i = 0; i < nAtoms; i++) {
        hash.withinRadii(xCoord[i], yCoord[i], zCoord[i], radius[i], neighbours)
        let ia = 0
        let ni = neighbours[ ia ]
        while (ni >= 0) {
            if (i < ni) {
            projectTorus(state, i, ni)
            }
            ni = neighbours[ ++ia ]
        }
    }
}

function fixNegatives (grid: NumberArray) {
    for (let i = 0; i < grid.length; i++) {
        if (grid[i] < 0) grid[i] = 0
    }
}

function fixAtomIDs (atomIndex: NumberArray, indexList: NumberArray) {
    for (let i = 0; i < atomIndex.length; i++) {
        atomIndex[i] = indexList[atomIndex[i]]
    }
}

//

interface MolSurfCalcState {
    /** Cached last value for obscured test */
    lastClip: number
    /** Neighbours as transient result array from lookup3d */
    neighbours: Result<number>

    lookup3D: Lookup3D
    position: PositionData
    radius: (index: number) => number
    delta: Vec3
    min: Vec3

    maxRadius: number

    n: number
    scaleFactor: number

    /** Angle lookup tables */
    cosTable: Float32Array
    sinTable: Float32Array

    probePositions: number
    ngTorus: number
}



export const MolecularSurfaceCalculationParams = {
    scaleFactor: PD.Numeric(2, { min: 0.1, max: 10, step: 0.1 }),
    probeRadius: PD.Numeric(1.4, { min: 0, max: 10, step: 0.1 }),
    probePositions: PD.Numeric(30, { min: 12, max: 90, step: 1 }),
}
export const DefaultMolecularSurfaceCalculationProps = PD.getDefaultValues(MolecularSurfaceCalculationParams)
export type MolecularSurfaceCalculationProps = typeof DefaultMolecularSurfaceCalculationProps

function createState(nAtoms: number, props: MolecularSurfaceCalculationProps): MolSurfCalcState {
    const { scaleFactor, probeRadius, probePositions } = props
    const { cosTable, sinTable } = getAngleTables(probePositions)
    const ngTorus = Math.max(5, 2 + Math.floor(probeRadius * scaleFactor))


    return {
        lastClip: -1,
        neighbours: Int32Array,

        xCoord: new Float32Array(nAtoms),
        yCoord: new Float32Array(nAtoms),
        zCoord: new Float32Array(nAtoms),
        radius: new Float32Array(nAtoms),
        radiusSq: new Float32Array(nAtoms),
        maxRadius: 0,

        n: nAtoms,
        scaleFactor,

        cosTable,
        sinTable,

        probePositions,
        ngTorus,
    }
}

//

export function MolecularSurface(coordList: Float32Array, radiusList: Float32Array, indexList: Uint16Array|Uint32Array) {
    // Field generation method adapted from AstexViewer (Mike Hartshorn)
    // by Fred Ludlow.
    // Other parts based heavily on NGL (Alexander Rose) EDT Surface class
    //
    // Should work as a drop-in alternative to EDTSurface (though some of
    // the EDT paramters are not relevant in this method).

    const nAtoms = radiusList.length

    const x = new Float32Array(nAtoms)
    const y = new Float32Array(nAtoms)
    const z = new Float32Array(nAtoms)

    for (let i = 0; i < nAtoms; i++) {
        const ci = 3 * i
        x[ i ] = coordList[ ci ]
        y[ i ] = coordList[ ci + 1 ]
        z[ i ] = coordList[ ci + 2 ]
    }

    let bbox = computeBoundingBox(coordList)
    if (coordList.length === 0) {
        bbox[ 0 ].set([ 0, 0, 0 ])
        bbox[ 1 ].set([ 0, 0, 0 ])
    }
    const min = bbox[0]
    const max = bbox[1]

    let r: Float32Array, r2: Float32Array // Atom positions, expanded radii (squared)
    let maxRadius: number

    // Parameters
    let probeRadius: number, scaleFactor: number, setAtomID: boolean, probePositions: number

    // Grid params
    let dim: Float32Array, matrix: Float32Array, grid: NumberArray, atomIndex: Int32Array

    // grid indices -> xyz coords
    let gridx: Float32Array, gridy: Float32Array, gridz: Float32Array

    // Spatial Hash
    let hash: iAVHash

    // Neighbour array to be filled by hash
    let neighbours: Int32Array

    let ngTorus: number

    function init (_probeRadius?: number, _scaleFactor?: number, _setAtomID?: boolean, _probePositions?: number) {
        probeRadius = defaults(_probeRadius, 1.4)
        scaleFactor = defaults(_scaleFactor, 2.0)
        setAtomID = defaults(_setAtomID, true)
        probePositions = defaults(_probePositions, 30)

        r = new Float32Array(nAtoms)
        r2 = new Float32Array(nAtoms)

        for (let i = 0; i < r.length; ++i) {
            var rExt = radiusList[ i ] + probeRadius
            r[ i ] = rExt
            r2[ i ] = rExt * rExt
        }

        maxRadius = 0
        for (let j = 0; j < r.length; ++j) {
            if (r[ j ] > maxRadius) maxRadius = r[ j ]
        }

        initializeGrid()
        getAngleTables(probePositions)
        initializeHash()

        lastClip = -1
    }

    function initializeGrid () {
        const surfGrid = getSurfaceGrid(
            min, max, maxRadius, scaleFactor, 0.0
        )

        scaleFactor = surfGrid.scaleFactor
        dim = surfGrid.dim
        matrix = surfGrid.matrix

        ngTorus = Math.max(5, 2 + Math.floor(probeRadius * scaleFactor))

        grid = fillUniform(new Float32Array(dim[0] * dim[1] * dim[2]), -1001.0)

        atomIndex = new Int32Array(grid.length)

        gridx = new Float32Array(dim[0])
        gridy = new Float32Array(dim[1])
        gridz = new Float32Array(dim[2])

        fillGridDim(gridx, min[0], 1 / scaleFactor)
        fillGridDim(gridy, min[1], 1 / scaleFactor)
        fillGridDim(gridz, min[2], 1 / scaleFactor)
    }



    function initializeHash () {
        hash = makeAVHash(x, y, z, r, min, max, 2.01 * maxRadius)
        neighbours = new Int32Array(hash.neighbourListLength)
    }





    function getVolume (probeRadius: number, scaleFactor: number, setAtomID: boolean) {
        // Basic steps are:
        // 1) Initialize
        // 2) Project points
        // 3) Project torii

        console.time('AVSurface.getVolume')

        console.time('AVSurface.init')
        init(probeRadius, scaleFactor, setAtomID)
        console.timeEnd('AVSurface.init')

        console.time('AVSurface.projectPoints')
        projectPoints()
        console.timeEnd('AVSurface.projectPoints')

        console.time('AVSurface.projectTorii')
        projectTorii()
        console.timeEnd('AVSurface.projectTorii')
        fixNegatives()
        fixAtomIDs()

        console.timeEnd('AVSurface.getVolume')
    }
}