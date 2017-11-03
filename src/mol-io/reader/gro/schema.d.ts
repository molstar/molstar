/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-base/collections/database'

export interface Header {
    title: string,
    timeInPs: number,
    /** number of decimal places */
    precision: { position: number, velocity: number },
    hasVelocities: boolean,
    box: [number, number, number]
}

export interface Atoms {
    count: number,
    residueNumber: Column<number>,
    residueName: Column<string>,
    atomName: Column<string>,
    atomNumber: Column<number>,
    x: Column<number>,
    y: Column<number>,
    z: Column<number>,
    vx: Column<number>,
    vy: Column<number>,
    vz: Column<number>
}

export interface Structure {
    header: Readonly<Header>,
    atoms: Readonly<Atoms>
}

export interface File {
    structures: Structure[]
}