/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from './color';

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
    colors: Color[]
}
export function ScaleLegend(minLabel: string, maxLabel: string, colors: Color[]): ScaleLegend {
    return { kind: 'scale-legend', minLabel, maxLabel, colors };
}