/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import * as React from 'react';
import { _Props, _State } from '../base';
import { Legend as LegendData, ScaleLegend as ScaleLegendData, TableLegend as TableLegendData } from '../../mol-util/legend';

export type LegendProps<L extends LegendData> = { legend: L }
export type Legend = React.ComponentClass<LegendProps<any>>

export function legendFor(legend: LegendData): Legend | undefined {
    switch (legend.kind) {
        case 'scale-legend': return ScaleLegend;
        case 'table-legend': return TableLegend;
        default:
            const _: never = legend;
            console.warn(`${_} has no associated UI component`);
            return void 0;
    }
}

export class ScaleLegend extends React.PureComponent<LegendProps<ScaleLegendData>> {
    render() {
        const { legend } = this.props;
        const colors = legend.colors.map(c => Color.toStyle(c)).join(', ');
        return  <div className='msp-scale-legend'>
            <div style={{ background: `linear-gradient(to right, ${colors})` }}>
                <span style={{float: 'left'}}>{legend.minLabel}</span>
                <span style={{float: 'right'}}>{legend.maxLabel}</span>
            </div>
        </div>;
    }
}

export class TableLegend extends React.PureComponent<LegendProps<TableLegendData>> {
    render() {
        const { legend } = this.props;
        return <div className='msp-table-legend'>
            {legend.table.map((value, i) => {
                const [name, color] = value;
                return <div key={i}>
                    <div className='msp-table-legend-color' style={{backgroundColor: Color.toStyle(color)}}></div>
                    <div className='msp-table-legend-text'>{name}</div>
                </div>;
            })}
        </div>;
    }
}