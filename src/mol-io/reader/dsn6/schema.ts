/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface Dsn6Header {
    /** the origin of the map in grid units, along X, Y, Z */
    xStart: number
    yStart: number
    zStart: number
    /** the extent (size) of the map, along X, Y, Z, in grid units */
    xExtent: number
    yExtent: number
    zExtent: number
    /** number of grid points along the whole unit cell, along X, Y, Z */
    xRate: number
    yRate: number
    zRate: number
    /** Unit cell parameters */
    xlen: number
    ylen: number
    zlen: number
    alpha: number
    beta: number
    gamma: number
    /**
     * Constants that bring the electron density from byte to normal scale.
     * They are calculated like this: prod = 255.0/(rhomax-rhomin), plus = -rhomin*prod.
     */
    divisor: number
    summand: number
    /** Rms deviation of electron density map (only given in BRIX but not in DSN6) */
    sigma: number | undefined
}

/**
 * DSN6 http://www.uoxray.uoregon.edu/tnt/manual/node104.html
 * BRIX http://svn.cgl.ucsf.edu/svn/chimera/trunk/libs/VolumeData/dsn6/brix-1.html
 */
export interface Dsn6File {
    name: string
    header: Dsn6Header
    values: Float32Array
}