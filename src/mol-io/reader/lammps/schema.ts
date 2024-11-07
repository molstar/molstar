/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Column } from '../../../mol-data/db';

interface LammpsUnitStyle {
    mass: string;
    distance: string;
    time: string;
    energy: string;
    velocity: string;
    force: string;
    torque: string;
    temperature: string;
    pressure: string;
    viscosity: string;
    charge: string;
    dipole?: string;
    electricField?: string;
    density: string;
    scale: number;
}

export const lammpsUnitStyles: { [key: string]: LammpsUnitStyle } = {
    lj: {
        mass: 'unitless',
        distance: 'unitless',
        time: 'unitless',
        energy: 'unitless',
        velocity: 'unitless',
        force: 'unitless',
        torque: 'unitless',
        temperature: 'unitless',
        pressure: 'unitless',
        viscosity: 'unitless',
        charge: 'unitless',
        density: 'unitless',
        scale: 1.0,
    },
    real: {
        mass: 'grams/mole',
        distance: 'Angstroms',
        time: 'femtoseconds',
        energy: 'Kcal/mol',
        velocity: 'Angstroms/femtosecond',
        force: 'Kcal/mol-Angstrom',
        torque: 'Kcal/mol',
        temperature: 'Kelvin',
        pressure: 'atmospheres',
        viscosity: 'Poise',
        charge: 'multiple of electron charge',
        dipole: 'charge*Angstroms',
        electricField: 'volts/Angstrom',
        density: 'g/cm^3',
        scale: 1.0,
    },
    metal: {
        mass: 'grams/mole',
        distance: 'Angstroms',
        time: 'picoseconds',
        energy: 'eV',
        velocity: 'Angstroms/picosecond',
        force: 'eV/Angstrom',
        torque: 'eV',
        temperature: 'Kelvin',
        pressure: 'bars',
        viscosity: 'Poise',
        charge: 'multiple of electron charge',
        dipole: 'charge*Angstroms',
        electricField: 'volts/Angstrom',
        density: 'g/cm^3',
        scale: 1.0,
    },
    si: {
        mass: 'kilograms',
        distance: 'meters',
        time: 'seconds',
        energy: 'Joules',
        velocity: 'meters/second',
        force: 'Newtons',
        torque: 'Newton-meters',
        temperature: 'Kelvin',
        pressure: 'Pascals',
        viscosity: 'Pascal*second',
        charge: 'Coulombs',
        dipole: 'Coulombs*meters',
        electricField: 'volts/meter',
        density: 'kg/m^3',
        scale: 1.0, // leave as is
    },
    cgs: {
        mass: 'grams',
        distance: 'centimeters',
        time: 'seconds',
        energy: 'ergs',
        velocity: 'centimeters/second',
        force: 'dynes',
        torque: 'dyne-centimeters',
        temperature: 'Kelvin',
        pressure: 'dyne/cm^2',
        viscosity: 'Poise',
        charge: 'statcoulombs',
        dipole: 'statcoul-cm',
        electricField: 'statvolt/cm',
        density: 'g/cm^3',
        scale: 1.0, // leave as is
    },
    electron: {
        mass: 'atomic mass units',
        distance: 'Bohr',
        time: 'femtoseconds',
        energy: 'Hartrees',
        velocity: 'Bohr/atomic time units',
        force: 'Hartrees/Bohr',
        temperature: 'Kelvin',
        pressure: 'Pascals',
        charge: 'multiple of electron charge',
        dipole: 'Debye',
        electricField: 'volts/cm',
        density: 'unitless',
        torque: '',
        viscosity: '',
        scale: 0.529177,
    },
    micro: {
        mass: 'picograms',
        distance: 'micrometers',
        time: 'microseconds',
        energy: 'picogram-micrometer^2/microsecond^2',
        velocity: 'micrometers/microsecond',
        force: 'picogram-micrometer/microsecond^2',
        torque: 'picogram-micrometer^2/microsecond^2',
        temperature: 'Kelvin',
        pressure: 'picogram/(micrometer-microsecond^2)',
        viscosity: 'picogram/(micrometer-microsecond)',
        charge: 'picocoulombs',
        dipole: 'picocoulomb-micrometer',
        electricField: 'volt/micrometer',
        density: 'pg/μm^3',
        scale: 1.0, // leave as is
    },
    nano: {
        mass: 'attograms',
        distance: 'nanometers',
        time: 'nanoseconds',
        energy: 'attogram-nanometer^2/nanosecond^2',
        velocity: 'nanometers/nanosecond',
        force: 'attogram-nanometer/nanosecond^2',
        torque: 'attogram-nanometer^2/nanosecond^2',
        temperature: 'Kelvin',
        pressure: 'attogram/(nanometer-nanosecond^2)',
        viscosity: 'attogram/(nanometer-nanosecond)',
        charge: 'multiple of electron charge',
        dipole: 'charge-nanometer',
        electricField: 'volt/nanometer',
        density: 'ag/nm^3',
        scale: 10.0,
    }
};

export const UnitStyles = ['real', 'metal', 'si', 'cgs', 'electron', 'micro', 'nano', 'lj'] as const;
export type UnitStyle = typeof UnitStyles[number];

export interface LammpsDataFile {
    readonly atoms: {
        readonly count: number
        readonly atomId: Column<number>
        readonly moleculeId: Column<number>
        readonly atomType: Column<number>
        readonly charge: Column<number>
        readonly x: Column<number>,
        readonly y: Column<number>,
        readonly z: Column<number>,
    }
    readonly bonds: {
        readonly count: number
        readonly bondId: Column<number>
        readonly bondType: Column<number>
        readonly atomIdA: Column<number>
        readonly atomIdB: Column<number>
    }
}

export interface LammpsBox {
    lower: [number, number, number],
    length: [number, number, number],
    periodicity: [string, string, string]
}

export interface LammpsFrame {
    count: number,
    atomMode: string,
    atomId: Column<number>,
    moleculeId: Column<number>,
    atomType: Column<number>,
    x: Column<number>,
    y: Column<number>,
    z: Column<number>,
}

export interface LammpsTrajectoryFile {
    frames: LammpsFrame[],
    times: number[],
    bounds: LammpsBox[],
    timeOffset: number,
    deltaTime: number
}
