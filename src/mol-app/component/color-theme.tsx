/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ColorTheme } from 'mol-theme/color';
import { Color } from 'mol-util/color';

export interface ColorThemeComponentProps {
    colorTheme: ColorTheme
}

export interface ColorThemeComponentState {

}

export class ColorThemeComponent extends React.Component<ColorThemeComponentProps, ColorThemeComponentState> {
    state = {

    }

    render() {
        const ct = this.props.colorTheme
        return <div>
            <span>Color Theme </span>

            {ct.description ? <div><i>{ct.description}</i></div> : ''}
            {
                ct.legend && ct.legend.kind === 'scale-legend'
                    ? <div
                        style={{
                            width: '100%',
                            height: '30px',
                            background: `linear-gradient(to right, ${ct.legend.colors.map(c => Color.toStyle(c)).join(', ')})`
                        }}
                    >
                        <span style={{float: 'left', padding: '6px', color: 'white', fontWeight: 'bold', backgroundColor: 'rgba(0, 0, 0, 0.2)'}}>{ct.legend.minLabel}</span>
                        <span style={{float: 'right', padding: '6px', color: 'white', fontWeight: 'bold', backgroundColor: 'rgba(0, 0, 0, 0.2)'}}>{ct.legend.maxLabel}</span>
                    </div>
                : ct.legend && ct.legend.kind === 'table-legend'
                    ? <div>
                        {ct.legend.table.map((value, i) => {
                            const [name, color] = value
                            return <div key={i} style={{minWidth: '60px', marginRight: '5px', display: 'inline-block'}}>
                                <div style={{width: '30px', height: '20px', backgroundColor: Color.toStyle(color), display: 'inline-block'}}></div>
                                {name}
                            </div>
                        })}
                    </div>
                : ''
            }
        </div>;
    }
}