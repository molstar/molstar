/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { DataLocation } from '../../../mol-model/location';
import { ElementIndex, Structure, StructureElement, Unit } from '../../../mol-model/structure';

export namespace ConfalPyramidsTypes {
    export type Pyramid = {
        PDB_model_number: number,
        name: string,
        auth_asym_id_1: string,
        auth_seq_id_1: number,
        label_comp_id_1: string,
        label_alt_id_1: string,
        PDB_ins_code_1: string,
        auth_asym_id_2: string,
        auth_seq_id_2: number,
        label_comp_id_2: string,
        label_alt_id_2: string,
        PDB_ins_code_2: string,
        confal_score: number,
        NtC: string
    }

    export interface PyramidsData {
        pyramids: Array<Pyramid>,
        names: Map<string, number>,
        locations: Array<Location>,
        hasMultipleModels: boolean
    }

    export interface LocationData {
        readonly pyramid: Pyramid
        readonly isLower: boolean;
    }

    export interface Element {
        structure: Structure;
        unit: Unit.Atomic;
        element: ElementIndex;
    }

    export interface Location extends DataLocation<LocationData, Element> {}

    export function Location(pyramid: Pyramid, isLower: boolean, structure?: Structure, unit?: Unit.Atomic, element?: ElementIndex) {
        return DataLocation('pyramid', { pyramid, isLower }, { structure: structure as any, unit: unit as any, element: element as any });
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'data-location' && x.tag === 'pyramid';
    }

    export function toElementLocation(location: Location) {
        return StructureElement.Location.create(location.element.structure, location.element.unit, location.element.element);
    }
}
