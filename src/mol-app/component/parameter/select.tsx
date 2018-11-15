/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from 'mol-util/param-definition';

export interface SelectParamComponentProps<T extends string> {
    label: string
    param: PD.Select<T>
    value: T
    onChange(v: T): void
}

export interface SelectParamComponentState<T extends string> {
    value: T
}

export class SelectParamComponent<T extends string> extends React.Component<SelectParamComponentProps<T>, SelectParamComponentState<T>> {
    state = {
        value: this.props.value
    }

    onChange(value: T) {
        this.setState({ value })
        this.props.onChange(value)
    }

    render() {
        return <div>
            <span>{this.props.label} </span>
            <select value={this.state.value} onChange={e => this.onChange(e.target.value as T) }>
                {this.props.param.options.map(v => {
                    const [value, label] = v
                    return <option key={label} value={value}>{label}</option>
                })}
            </select>
        </div>;
    }
}