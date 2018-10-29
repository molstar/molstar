/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { BooleanParam } from 'mol-util/parameter';

export interface BooleanParamComponentProps {
    param: BooleanParam
    value: boolean
    onChange(v: boolean): void
}

export interface BooleanParamComponentState {
    value: boolean
}

export class BooleanParamComponent extends React.Component<BooleanParamComponentProps, BooleanParamComponentState> {
    state = {
        value: this.props.value
    }

    onChange(value: boolean) {
        this.setState({ value })
        this.props.onChange(value)
    }

    render() {
        return <div>
            <span>{this.props.param.label} </span>
            <button onClick={e => this.onChange(!this.state.value) }>
                {this.state.value ? 'Off' : 'On'}
            </button>
        </div>;
    }
}