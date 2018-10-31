/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { MultiSelectParam } from 'mol-util/parameter';

export interface MultiSelectParamComponentProps<T extends string> {
    param: MultiSelectParam<T>
    value: T[]
    onChange(v: T[]): void
}

export interface MultiSelectParamComponentState<T extends string> {
    value: T[]
}

export class MultiSelectParamComponent<T extends string> extends React.Component<MultiSelectParamComponentProps<T>, MultiSelectParamComponentState<T>> {
    state = {
        value: this.props.value
    }

    onChange(value: T[]) {
        this.setState({ value })
        this.props.onChange(value)
    }

    render() {
        return <div>
            <span>{this.props.param.label} </span>
            <select multiple value={this.state.value} onChange={e => {
                const value = Array.from(e.target.options).filter(option => option.selected).map(option => option.value)
                this.onChange(value as T[]) 
            }}>
                {this.props.param.options.map(v => {
                    const [value, label] = v
                    return <option key={label} value={value}>{label}</option>
                })}
            </select>
        </div>;
    }
}