/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { TextParam } from 'mol-view/parameter';

export interface TextParamComponentProps {
    param: TextParam
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
            <span>{this.props.param.label} </span>
            <input type='text'
                value={this.state.value}
                onChange={e => this.onChange(e.currentTarget.value)}
            >
            </input>
        </div>;
    }
}