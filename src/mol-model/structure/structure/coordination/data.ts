/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from '../unit';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ElementIndex } from '../../model/indexing';
import { StructureElement } from '../element';

export interface CoordinationSite {
    /** Central atom's unit */
    readonly unit: Unit.Atomic
    /** Central atom's element index */
    readonly element: ElementIndex
    /** Central atom's UnitIndex within its unit */
    readonly unitIndex: StructureElement.UnitIndex
    /** Positions of coordinating neighbor atoms */
    readonly ligandPositions: ReadonlyArray<Vec3>
    /** Neighbor atom units */
    readonly ligandUnits: ReadonlyArray<Unit.Atomic>
    /** Neighbor atom element indices */
    readonly ligandElements: ReadonlyArray<ElementIndex>
}

export interface Coordination {
    readonly sites: ReadonlyArray<CoordinationSite>
    /** Get site indices for atoms in the given unit */
    readonly getSiteIndices: (unit: Unit.Atomic, element: ElementIndex) => ReadonlyArray<number>
}

const EmptyArray: ReadonlyArray<any> = [];

export const EmptyCoordination: Coordination = {
    sites: EmptyArray,
    getSiteIndices: () => EmptyArray,
};
