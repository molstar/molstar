/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { DnatcoTypes } from '../types';
import { DataLocation } from '../../../mol-model/location';
import { DataLoci } from '../../../mol-model/loci';
import { confalPyramidLabel } from './behavior';

export namespace ConfalPyramidsTypes {
    export interface Location extends DataLocation<DnatcoTypes.HalfStep, {}> {}

    export function Location(step: DnatcoTypes.Step, isLower: boolean) {
        return DataLocation(DnatcoTypes.DataTag, { step, isLower }, {});
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'data-location' && x.tag === DnatcoTypes.DataTag;
    }

    export interface Loci extends DataLoci<DnatcoTypes.Step[], number> {}

    export function Loci(data: DnatcoTypes.Step[], elements: ReadonlyArray<number>): Loci {
        return DataLoci(DnatcoTypes.DataTag, data, elements, undefined, () => elements[0] !== undefined ? confalPyramidLabel(data[elements[0]]) : '');
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'data-loci' && x.tag === DnatcoTypes.DataTag;
    }
}
