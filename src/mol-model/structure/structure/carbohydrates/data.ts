/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../unit';
import { Vec3 } from 'mol-math/linear-algebra';
import { ResidueIndex, ElementIndex } from '../../model';
import { SaccharideComponent } from './constants';
import StructureElement from '../element';

export interface CarbohydrateLink {
    readonly carbohydrateIndexA: number
    readonly carbohydrateIndexB: number
}

export interface CarbohydrateTerminalLink {
    readonly carbohydrateIndex: number
    readonly elementIndex: StructureElement.UnitIndex
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

/** partial carbohydrate with no ring present */
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
    getLinkIndices: (unit: Unit, anomericCarbon: ElementIndex) => ReadonlyArray<number>
    getTerminalLinkIndex: (unitA: Unit, elementA: ElementIndex, unitB: Unit, elementB: ElementIndex) => number | undefined
    getTerminalLinkIndices: (unit: Unit, element: ElementIndex) => ReadonlyArray<number>
    getAnomericCarbons: (unit: Unit, residueIndex: ResidueIndex) => ReadonlyArray<ElementIndex>
}

const EmptyArray: ReadonlyArray<any> = []
export const EmptyCarbohydrates: Carbohydrates = {
    links: EmptyArray,
    terminalLinks: EmptyArray,
    elements: EmptyArray,
    partialElements: EmptyArray,
    getElementIndex: () => undefined,
    getLinkIndex: () => undefined,
    getLinkIndices: () => EmptyArray,
    getTerminalLinkIndex: () => undefined,
    getTerminalLinkIndices: () => EmptyArray,
    getAnomericCarbons: () => [],
}