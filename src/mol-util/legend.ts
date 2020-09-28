/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from './color';
import { ColorListEntry } from './color/color';

export type Legend = TableLegend | ScaleLegend

export interface TableLegend {
    kind: 'table-legend'
    table: [ string, Color ][]
}
export function TableLegend(table: [ string, Color ][]): TableLegend {
    return { kind: 'table-legend', table };
}

export interface ScaleLegend {
    kind: 'scale-legend'
    minLabel: string,
    maxLabel: string,
    colors: ColorListEntry[]
}
export function ScaleLegend(minLabel: string, maxLabel: string, colors: ColorListEntry[]): ScaleLegend {
    return { kind: 'scale-legend', minLabel, maxLabel, colors };
}