/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { BooleanParamComponent } from './parameter/boolean';
import { NumberParamComponent } from './parameter/number';
import { RangeParamComponent } from './parameter/range';
import { SelectParamComponent } from './parameter/select';
import { MultiSelectParamComponent } from './parameter/multi-select';
import { TextParamComponent } from './parameter/text';
import { ColorParamComponent } from './parameter/color';

interface ParametersProps<P extends PD.Params> {
    params: P
    values: { [k in keyof P]: P[k]['defaultValue'] }
    onChange<K extends keyof P>(k: K, v: P[K]['defaultValue']): void
}

type ParametersState = {}

function getParamComponent<P extends PD.Any>(p: PD.Any, value: P['defaultValue'], onChange: (v: P['defaultValue']) => void) {
    switch (p.type) {
        case 'boolean':
            return <BooleanParamComponent param={p} value={value} onChange={onChange} />
        case 'number':
            return <NumberParamComponent param={p} value={value} onChange={onChange} />
        case 'range':
            return <RangeParamComponent param={p} value={value} onChange={onChange} />
        case 'select':
            return <SelectParamComponent param={p} value={value} onChange={onChange} />
        case 'multi-select':
            return <MultiSelectParamComponent param={p} value={value} onChange={onChange} />
        case 'text':
            return <TextParamComponent param={p} value={value} onChange={onChange} />
        case 'color':
            return <ColorParamComponent param={p} value={value} onChange={onChange} />
    }
    return ''
}

export class ParametersComponent<P extends PD.Params> extends React.Component<ParametersProps<P>, ParametersState> {
    onChange(k: string, value: any) {
        this.props.onChange(k, value)
    }

    render() {
        return <div style={{ width: '100%' }}>
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