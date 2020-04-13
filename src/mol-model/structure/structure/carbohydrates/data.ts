/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../unit';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ResidueIndex, ElementIndex } from '../../model';
import { SaccharideComponent } from './constants';
import StructureElement from '../element';
import { UnitRings } from '../unit/rings';

export interface CarbohydrateLink {
    readonly carbohydrateIndexA: number
    readonly carbohydrateIndexB: number
}

export interface CarbohydrateTerminalLink {
    readonly carbohydrateIndex: number
    readonly elementIndex: StructureElement.UnitIndex
    readonly elementUnit: Unit.Atomic
    /** specifies direction of the link */
    readonly fromCarbohydrate: boolean
}

export interface CarbohydrateElement {
    readonly geometry: { readonly center: Vec3, readonly normal: Vec3, readonly direction: Vec3 },
    readonly unit: Unit.Atomic,
    readonly residueIndex: ResidueIndex,
    readonly component: SaccharideComponent,
    readonly ringIndex: UnitRings.Index,
    readonly altId: string
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

    getElementIndices: (unit: Unit.Atomic, element: ElementIndex) => ReadonlyArray<number>
    getLinkIndices: (unit: Unit.Atomic, element: ElementIndex) => ReadonlyArray<number>
    getTerminalLinkIndices: (unit: Unit.Atomic, element: ElementIndex) => ReadonlyArray<number>
}

const EmptyArray: ReadonlyArray<any> = [];
export const EmptyCarbohydrates: Carbohydrates = {
    links: EmptyArray,
    terminalLinks: EmptyArray,
    elements: EmptyArray,
    partialElements: EmptyArray,

    getElementIndices: () => EmptyArray,
    getLinkIndices: () => EmptyArray,
    getTerminalLinkIndices: () => EmptyArray,
};