/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { BooleanParamComponent } from './parameter/boolean';
import { NumberParamComponent } from './parameter/number';
import { SelectParamComponent } from './parameter/select';
import { MultiSelectParamComponent } from './parameter/multi-select';
import { TextParamComponent } from './parameter/text';
import { ColorParamComponent } from './parameter/color';
import { camelCaseToWords } from 'mol-util/string';

interface ParametersProps<P extends PD.Params> {
    params: P
    values: { [k in keyof P]: P[k]['defaultValue'] }
    onChange<K extends keyof P>(k: K, v: P[K]['defaultValue']): void
}

type ParametersState = {}

function getParamComponent<P extends PD.Any>(label: string, param: PD.Any, value: P['defaultValue'], onChange: (v: P['defaultValue']) => void) {
    switch (param.type) {
        case 'boolean':
            return <BooleanParamComponent label={label} param={param} value={value} onChange={onChange} />
        case 'number':
            return <NumberParamComponent label={label} param={param} value={value} onChange={onChange} />
        case 'select':
            return <SelectParamComponent label={label} param={param} value={value} onChange={onChange} />
        case 'multi-select':
            return <MultiSelectParamComponent label={label} param={param} value={value} onChange={onChange} />
        case 'text':
            return <TextParamComponent label={label} param={param} value={value} onChange={onChange} />
        case 'color':
            return <ColorParamComponent label={label} param={param} value={value} onChange={onChange} />
    }
    return ''
}

function getLabel(name: string, param: PD.Base<any>) {
    return param.label === undefined ? camelCaseToWords(name) : param.label
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
                const label = getLabel(k, param)
                return <div key={k}>
                    {getParamComponent(label, param, value, v => this.onChange(k, v))}
                </div>
            })}
        </div>;
    }
}