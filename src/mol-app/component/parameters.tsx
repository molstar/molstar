/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { Param, Params } from 'mol-view/parameter';
import { BooleanParamComponent } from './parameter/boolean';
import { NumberParamComponent } from './parameter/number';
import { RangeParamComponent } from './parameter/range';
import { SelectParamComponent } from './parameter/select';

interface ParametersProps<P extends Params> {
    params: P
    values: { [k in keyof P]: P[k]['defaultValue'] }
    onChange<K extends keyof P>(k: K, v: P[K]['defaultValue']): void
}

type ParametersState = {}

function getParamComponent<P extends Param>(p: Param, value: P['defaultValue'], onChange: (v: P['defaultValue']) => void) {
    switch (p.type) {
        case 'boolean':
            return <BooleanParamComponent param={p} value={value} onChange={onChange} />
        case 'number':
            return <NumberParamComponent param={p} value={value} onChange={onChange} />
        case 'range':
            return <RangeParamComponent param={p} value={value} onChange={onChange} />
        case 'select':
            return <SelectParamComponent param={p} value={value} onChange={onChange} />
    }
    return ''
}

export class ParametersComponent<P extends Params> extends React.Component<ParametersProps<P>, ParametersState> {
    onChange(k: string, value: any) {
        this.props.onChange(k, value)
    }

    render() {
        return <div>
            { Object.keys(this.props.params).map(k => {
                const param = this.props.params[k]
                const value = this.props.values[k]
                return <div key={k}>
                    {getParamComponent(param, value, v => this.onChange(k, v))}
                </div>
            })}
        </div>;
    }
}