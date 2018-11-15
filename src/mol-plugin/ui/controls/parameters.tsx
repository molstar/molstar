/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { camelCaseToWords } from 'mol-util/string';

export interface ParameterControlsProps<P extends PD.Params = PD.Params> {
    params: P,
    values: any,
    onChange: ParamOnChange,
    isEnabled?: boolean,
    onEnter?: () => void
}

export class ParameterControls<P extends PD.Params> extends React.PureComponent<ParameterControlsProps<P>, {}> {
    render() {
        const common = {
            onChange: this.props.onChange,
            isEnabled: this.props.isEnabled,
            onEnter: this.props.onEnter,
        }
        const params = this.props.params;
        const values = this.props.values;
        return <div style={{ width: '100%' }}>
            {Object.keys(params).map(key => {
                const param = params[key];
                if (param.type === 'value') return null;
                if (param.type === 'mapped') return <MappedControl param={param} key={key} {...common} name={key} value={values[key]} />
                if (param.type === 'group') return <GroupControl param={param} key={key} {...common} name={key} value={values[key]} />
                return <ParamWrapper control={controlFor(param)} param={param} key={key} {...common} name={key} value={values[key]} />
            })}
        </div>;
    }
}

function controlFor(param: PD.Any): ValueControl {
    switch (param.type) {
        case 'boolean': return BoolControl;
        case 'number': return NumberControl;
        case 'multi-select': return MultiSelectControl;
        case 'color': return ColorControl;
        case 'select': return SelectControl;
        case 'text': return TextControl;
        case 'interval': return IntervalControl;
        case 'converted': return ConvertedControl;
        case 'group': throw Error('Must be handled separately');
        case 'mapped': throw Error('Must be handled separately');
    }
    throw new Error('not supported');
}

type ParamWrapperProps = { name: string, value: any, param: PD.Base<any>, onChange: ParamOnChange, control: ValueControl, onEnter?: () => void, isEnabled?: boolean }
export type ParamOnChange = (params: { param: PD.Base<any>, name: string, value: any }) => void
type ValueControlProps<P extends PD.Base<any> = PD.Base<any>> = { value: any, param: P, isEnabled?: boolean, onChange: (v: any) => void, onEnter?: () => void }
type ValueControl = React.ComponentClass<ValueControlProps<any>>

function getLabel(name: string, param: PD.Base<any>) {
    return param.label === undefined ? camelCaseToWords(name) : param.label
}

export class ParamWrapper extends React.PureComponent<ParamWrapperProps> {
    onChange = (value: any) => {
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    render() {
        return <div style={{ padding: '0 3px', borderBottom: '1px solid #ccc' }}>
            <div style={{ lineHeight: '20px', float: 'left' }} title={this.props.param.description}>{getLabel(this.props.name, this.props.param)}</div>
            <div style={{ float: 'left', marginLeft: '5px' }}>
                <this.props.control value={this.props.value} param={this.props.param} onChange={this.onChange} onEnter={this.props.onEnter} isEnabled={this.props.isEnabled} />
            </div>
            <div style={{ clear: 'both' }} />
        </div>;
    }
}

export class BoolControl extends React.PureComponent<ValueControlProps> {
    onClick = () => {
        this.props.onChange(!this.props.value);
    }

    render() {
        return <button onClick={this.onClick} disabled={!this.props.isEnabled}>{this.props.value ? '✓ On' : '✗ Off'}</button>;
    }
}

export class NumberControl extends React.PureComponent<ValueControlProps<PD.Numeric>, { value: string }> {
    // state = { value: this.props.value }
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        this.props.onChange(+e.target.value);
        // this.setState({ value: e.target.value });
    }

    render() {
        return <input type='range'
            value={'' + this.props.value}
            min={this.props.param.min}
            max={this.props.param.max}
            step={this.props.param.step}
            onChange={this.onChange}
        />;
    }
}

export class TextControl extends React.PureComponent<ValueControlProps<PD.Text>> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const value = e.target.value;
        if (value !== this.props.value) {
            this.props.onChange(value);
        }
    }

    onKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if (!this.props.onEnter) return;
        if ((e.keyCode === 13 || e.charCode === 13)) {
            this.props.onEnter();
        }
    }

    render() {
        return <input type='text'
            value={this.props.value || ''}
            onChange={this.onChange}
            onKeyPress={this.props.onEnter ? this.onKeyPress : void 0}
        />;
    }
}

export class SelectControl extends React.PureComponent<ValueControlProps<PD.Select<any>>> {
    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        this.setState({ value: e.target.value });
        this.props.onChange(e.target.value);
    }

    render() {
        return <select value={this.props.value || ''} onChange={this.onChange}>
            {this.props.param.options.map(([value, label]) => <option key={value} value={value}>{label}</option>)}
        </select>;
    }
}

export class MultiSelectControl extends React.PureComponent<ValueControlProps<PD.MultiSelect<any>>> {
    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        const value = Array.from(e.target.options).filter(option => option.selected).map(option => option.value);
        this.setState({ value });
        this.props.onChange(value);
    }

    render() {
        return <select multiple value={this.props.value || ''} onChange={this.onChange}>
            {this.props.param.options.map(([value, label]) => <option key={label} value={value}>{label}</option>)}
        </select>;
    }
}

export class IntervalControl extends React.PureComponent<ValueControlProps<PD.Interval>> {
    // onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    //     this.setState({ value: e.target.value });
    //     this.props.onChange(e.target.value);
    // }

    render() {
        return <span>interval TODO</span>;
    }
}

export class ColorControl extends React.PureComponent<ValueControlProps<PD.Color>> {
    // onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    //     this.setState({ value: e.target.value });
    //     this.props.onChange(e.target.value);
    // }

    render() {
        return <span>color TODO</span>;
    }
}

export class ConvertedControl extends React.PureComponent<ValueControlProps<PD.Converted<any, any>>> {
    onChange = (v: any) => {
        console.log('onChange', v)
        this.props.onChange(this.props.param.toValue(v));
    }

    render() {
        const Control: ValueControl = controlFor(this.props.param.param as PD.Any);

        return <Control value={this.props.param.fromValue(this.props.value)} param={this.props.param.param} onChange={this.onChange} onEnter={this.props.onEnter} isEnabled={this.props.isEnabled} />
    }
}

type GroupWrapperProps = { name: string, value: PD.Group<any>['defaultValue'], param: PD.Group<any>, onChange: ParamOnChange, onEnter?: () => void, isEnabled?: boolean }
export class GroupControl extends React.PureComponent<GroupWrapperProps> {
    change(value: PD.Mapped<any>['defaultValue'] ) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeParam: ParamOnChange = e => {
        const value: PD.Mapped<any>['defaultValue'] = this.props.value;
        this.change({ ...value.params, [e.name]: e.value });
    }

    render() {
        const value: PD.Mapped<any>['defaultValue'] = this.props.value;
        const params = this.props.param.params;

        return <div>
            <ParameterControls params={params} onChange={this.onChangeParam} values={value.params} onEnter={this.props.onEnter} isEnabled={this.props.isEnabled} />
        </div>
    }
}

type MappedWrapperProps = { name: string, value: PD.Mapped<any>['defaultValue'], param: PD.Mapped<any>, onChange: ParamOnChange, onEnter?: () => void, isEnabled?: boolean }
export class MappedControl extends React.PureComponent<MappedWrapperProps> {
    change(value: PD.Mapped<any>['defaultValue'] ) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeName: ParamOnChange = e => {
        this.change({ name: e.value, params: this.props.param.map(e.value).defaultValue });
    }

    onChangeParam: ParamOnChange = e => {
        const value: PD.Mapped<any>['defaultValue'] = this.props.value;
        this.change({ name: value.name, params: e.value });
    }

    render() {
        const value: PD.Mapped<any>['defaultValue'] = this.props.value;
        const param = this.props.param.map(value.name);

        return <div>
            <ParamWrapper control={SelectControl} param={this.props.param.select}
                isEnabled={this.props.isEnabled} onChange={this.onChangeName} onEnter={this.props.onEnter}
                name={'name'} value={value.name} />
            <div style={{ borderLeft: '5px solid #777', paddingLeft: '5px' }}>
                {param.type === 'group'
                ? <GroupControl param={param} value={value} name='param' onChange={this.onChangeParam} onEnter={this.props.onEnter} isEnabled={this.props.isEnabled} />
                : param.type === 'mapped' || param.type === 'value' ? null
                : <ParamWrapper control={controlFor(param)} param={param} onChange={this.onChangeParam} onEnter={this.props.onEnter} isEnabled={this.props.isEnabled} name={'value'} value={value} />}
            </div>
        </div>
    }
}