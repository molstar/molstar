/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/tables';
import { Color } from '../../../mol-util/color';

export interface ColorParamComponentProps {
    label: string
    param: PD.Color
    value: Color
    onChange(v: Color): void
}

export interface ColorParamComponentState {
    value: Color
}

export class ColorParamComponent extends React.Component<ColorParamComponentProps, ColorParamComponentState> {
    state = {
        value: this.props.value
    }

    onChange(value: Color) {
        this.setState({ value })
        this.props.onChange(value)
    }

    render() {
        return <div>
            <span>{this.props.label} </span>
            <select value={this.state.value} onChange={e => this.onChange(Color(parseInt(e.target.value))) }>
                {Object.keys(ColorNames).map(name => {
                    return <option key={name} value={(ColorNames as { [k: string]: Color})[name]}>{name}</option>
                })}
            </select>
        </div>;
    }
}