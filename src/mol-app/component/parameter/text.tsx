/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from 'mol-util/param-definition';

export interface TextParamComponentProps {
    label: string
    param: PD.Text
    value: string
    onChange(v: string): void
}

export interface TextParamComponentState {
    value: string
}

export class TextParamComponent extends React.Component<TextParamComponentProps, TextParamComponentState> {
    state = {
        value: this.props.value
    }

    onChange(value: string) {
        this.setState({ value })
        this.props.onChange(value)
    }

    render() {
        return <div>
            <span>{this.props.label} </span>
            <input type='text'
                value={this.state.value}
                onChange={e => this.onChange(e.currentTarget.value)}
            >
            </input>
        </div>;
    }
}