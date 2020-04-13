/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../mol-data/db';

export interface GroHeader {
    title: string,
    timeInPs: number,
    /** number of decimal places */
    precision: { position: number, velocity: number },
    hasVelocities: boolean,
    box: [number, number, number]
}

export interface GroAtoms {
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

export interface GroStructure {
    header: Readonly<GroHeader>,
    atoms: Readonly<GroAtoms>
}

export interface GroFile {
    structures: GroStructure[]
}