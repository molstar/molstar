/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { DnatcoCommon as DC } from '../common';
import { DataLocation } from '../../../mol-model/location';
import { ElementIndex, Structure, StructureElement, Unit } from '../../../mol-model/structure';

export namespace ConfalPyramidsTypes {
    export type ConformerClasses = 'A' | 'B' | 'BII' | 'miB' | 'Z' | 'IC' | 'OPN' | 'SYN' | 'N';

    export const NtCToClasses: ReadonlyMap<string, [ConformerClasses, ConformerClasses]> = new Map([
        ['NANT', ['N', 'N']],
        ['AA00', ['A', 'A']],
        ['AA02', ['A', 'A']],
        ['AA03', ['A', 'A']],
        ['AA04', ['A', 'A']],
        ['AA08', ['A', 'A']],
        ['AA09', ['A', 'A']],
        ['AA01', ['A', 'A']],
        ['AA05', ['A', 'A']],
        ['AA06', ['A', 'A']],
        ['AA10', ['A', 'A']],
        ['AA11', ['A', 'A']],
        ['AA07', ['A', 'A']],
        ['AA12', ['A', 'A']],
        ['AA13', ['A', 'A']],
        ['AB01', ['A', 'B']],
        ['AB02', ['A', 'B']],
        ['AB03', ['A', 'B']],
        ['AB04', ['A', 'B']],
        ['AB05', ['A', 'B']],
        ['BA01', ['B', 'A']],
        ['BA05', ['B', 'A']],
        ['BA09', ['B', 'A']],
        ['BA08', ['BII', 'A']],
        ['BA10', ['B', 'A']],
        ['BA13', ['BII', 'A']],
        ['BA16', ['BII', 'A']],
        ['BA17', ['BII', 'A']],
        ['BB00', ['B', 'B']],
        ['BB01', ['B', 'B']],
        ['BB17', ['B', 'B']],
        ['BB02', ['B', 'B']],
        ['BB03', ['B', 'B']],
        ['BB11', ['B', 'B']],
        ['BB16', ['B', 'B']],
        ['BB04', ['B', 'BII']],
        ['BB05', ['B', 'BII']],
        ['BB07', ['BII', 'BII']],
        ['BB08', ['BII', 'BII']],
        ['BB10', ['miB', 'miB']],
        ['BB12', ['miB', 'miB']],
        ['BB13', ['miB', 'miB']],
        ['BB14', ['miB', 'miB']],
        ['BB15', ['miB', 'miB']],
        ['BB20', ['miB', 'miB']],
        ['IC01', ['IC', 'IC']],
        ['IC02', ['IC', 'IC']],
        ['IC03', ['IC', 'IC']],
        ['IC04', ['IC', 'IC']],
        ['IC05', ['IC', 'IC']],
        ['IC06', ['IC', 'IC']],
        ['IC07', ['IC', 'IC']],
        ['OP01', ['OPN', 'OPN']],
        ['OP02', ['OPN', 'OPN']],
        ['OP03', ['OPN', 'OPN']],
        ['OP04', ['OPN', 'OPN']],
        ['OP05', ['OPN', 'OPN']],
        ['OP06', ['OPN', 'OPN']],
        ['OP07', ['OPN', 'OPN']],
        ['OP08', ['OPN', 'OPN']],
        ['OP09', ['OPN', 'OPN']],
        ['OP10', ['OPN', 'OPN']],
        ['OP11', ['OPN', 'OPN']],
        ['OP12', ['OPN', 'OPN']],
        ['OP13', ['OPN', 'OPN']],
        ['OP14', ['OPN', 'OPN']],
        ['OP15', ['OPN', 'OPN']],
        ['OP16', ['OPN', 'OPN']],
        ['OP17', ['OPN', 'OPN']],
        ['OP18', ['OPN', 'OPN']],
        ['OP19', ['OPN', 'OPN']],
        ['OP20', ['OPN', 'OPN']],
        ['OP21', ['OPN', 'OPN']],
        ['OP22', ['OPN', 'OPN']],
        ['OP23', ['OPN', 'OPN']],
        ['OP24', ['OPN', 'OPN']],
        ['OP25', ['OPN', 'OPN']],
        ['OP26', ['OPN', 'OPN']],
        ['OP27', ['OPN', 'OPN']],
        ['OP28', ['OPN', 'OPN']],
        ['OP29', ['OPN', 'OPN']],
        ['OP30', ['OPN', 'OPN']],
        ['OP31', ['OPN', 'OPN']],
        ['OPS1', ['OPN', 'OPN']],
        ['OP1S', ['OPN', 'SYN']],
        ['AAS1', ['SYN', 'A']],
        ['AB1S', ['A', 'SYN']],
        ['AB2S', ['A', 'SYN']],
        ['BB1S', ['B', 'SYN']],
        ['BB2S', ['B', 'SYN']],
        ['BBS1', ['SYN', 'B']],
        ['ZZ01', ['Z', 'Z']],
        ['ZZ02', ['Z', 'Z']],
        ['ZZ1S', ['Z', 'SYN']],
        ['ZZ2S', ['Z', 'SYN']],
        ['ZZS1', ['SYN', 'Z']],
        ['ZZS2', ['SYN', 'Z']],
    ]);

    export type Pyramid = DC.NtCObject;

    export interface PyramidsData {
        pyramids: Array<Pyramid>;
        names: Map<string, number>;
        locations: Array<Location>;
        hasMultipleModels: boolean;
    }

    export interface LocationData {
        readonly pyramid: Pyramid;
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
