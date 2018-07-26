/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../unit';
import { Vec3 } from 'mol-math/linear-algebra';
import { ResidueIndex } from '../../model';
import { SaccharideComponent } from './constants';

export interface CarbohydrateLink {
    readonly carbohydrateIndexA: number
    readonly carbohydrateIndexB: number
}

export interface CarbohydrateTerminalLink {
    readonly carbohydrateIndex: number
    readonly elementIndex: number
    readonly elementUnit: Unit
    /** specifies direction of the link */
    readonly fromCarbohydrate: boolean
}

export interface CarbohydrateElement {
    readonly center: Vec3,
    readonly normal: Vec3,
    readonly direction: Vec3,
    readonly unit: Unit.Atomic
    readonly residueIndex: ResidueIndex
    readonly component: SaccharideComponent
}

export interface Carbohydrates {
    links: ReadonlyArray<CarbohydrateLink>
    terminalLinks: ReadonlyArray<CarbohydrateTerminalLink>
    elements: ReadonlyArray<CarbohydrateElement>
}