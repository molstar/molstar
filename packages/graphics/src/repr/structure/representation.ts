/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import type { StructureUnitTransforms } from '@molstar/model/model/structure/structure/util/unit-transforms';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Representation, type RepresentationProps, type RepresentationProvider } from '../representation.js';

export interface StructureRepresentationState extends Representation.State {
    unitTransforms: StructureUnitTransforms | null,
    unitTransformsVersion: number
}

export const StructureRepresentationStateBuilder: Representation.StateBuilder<StructureRepresentationState> = {
    create: () => {
        return {
            ...Representation.createState(),
            unitTransforms: null,
            unitTransformsVersion: -1
        };
    },
    update: (state: StructureRepresentationState, update: Partial<StructureRepresentationState>) => {
        Representation.updateState(state, update);
        if (update.unitTransforms !== undefined) state.unitTransforms = update.unitTransforms;
    }
};

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P, StructureRepresentationState> { }

export type StructureRepresentationProvider<P extends PD.Params, Id extends string = string> = RepresentationProvider<Structure, P, StructureRepresentationState, Id>
export function StructureRepresentationProvider<P extends PD.Params, Id extends string>(p: StructureRepresentationProvider<P, Id>): StructureRepresentationProvider<P, Id> { return p; }

//

export { ComplexRepresentation } from './complex-representation.js';
export { ComplexVisual } from './complex-visual.js';
export { UnitsRepresentation } from './units-representation.js';
export { UnitsVisual } from './units-visual.js';
