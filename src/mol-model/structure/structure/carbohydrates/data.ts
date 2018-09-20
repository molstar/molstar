/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../unit';
import { Vec3 } from 'mol-math/linear-algebra';
import { ResidueIndex, ElementIndex } from '../../model';
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
    readonly geometry: { readonly center: Vec3, readonly normal: Vec3, readonly direction: Vec3 },
    readonly anomericCarbon: ElementIndex,
    readonly unit: Unit.Atomic,
    readonly residueIndex: ResidueIndex,
    readonly component: SaccharideComponent,
    readonly ringAltId: string,
}

// partial carbohydrate with no ring present
export interface PartialCarbohydrateElement {
    readonly unit: Unit.Atomic,
    readonly residueIndex: ResidueIndex,
    readonly component: SaccharideComponent,
}

export interface Carbohydrates {
    links: ReadonlyArray<CarbohydrateLink>
    terminalLinks: ReadonlyArray<CarbohydrateTerminalLink>
    elements: ReadonlyArray<CarbohydrateElement>
    partialElements: ReadonlyArray<PartialCarbohydrateElement>
    getElementIndex: (unit: Unit, anomericCarbon: ElementIndex) => number | undefined
    getLinkIndex: (unitA: Unit, anomericCarbonA: ElementIndex, unitB: Unit, anomericCarbonB: ElementIndex) => number | undefined
    getAnomericCarbon: (unit: Unit, residueIndex: ResidueIndex) => ElementIndex | undefined
}