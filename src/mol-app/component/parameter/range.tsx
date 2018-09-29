/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { RangeParam } from 'mol-view/parameter';

export interface RangeParamComponentProps {
    param: RangeParam
    value: number
    onChange(v: number): void
}

export interface RangeParamComponentState {
    value: number
}

export class RangeParamComponent extends React.Component<RangeParamComponentProps, RangeParamComponentState> {
    state = {
        value: this.props.value
    }

    onChange(valueStr: string) {
        const value = Number.isInteger(this.props.param.step) ? parseInt(valueStr) : parseFloat(valueStr)
        this.setState({ value })
        this.props.onChange(value)
    }

    render() {
        return <div>
            <span>{this.props.param.label} </span>
            <input type='range'
                value={this.state.value}
                min={this.props.param.min}
                max={this.props.param.max}
                step={this.props.param.step}
                onChange={e => this.onChange(e.currentTarget.value)}
            >
            </input>
        </div>;
    }
}