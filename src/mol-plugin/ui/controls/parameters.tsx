/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { camelCaseToWords } from 'mol-util/string';
import { ColorNames } from 'mol-util/color/tables';
import { Color } from 'mol-util/color';

export interface ParameterControlsProps<P extends PD.Params = PD.Params> {
    params: P,
    values: any,
    onChange: ParamOnChange,
    isDisabled?: boolean,
    onEnter?: () => void
}

export class ParameterControls<P extends PD.Params> extends React.PureComponent<ParameterControlsProps<P>, {}> {
    render() {
        const params = this.props.params;
        const values = this.props.values;
        return <div style={{ width: '100%' }}>
            {Object.keys(params).map(key => {
                const param = params[key];
                const Control = controlFor(param);
                if (!Control) return null;
                return <Control param={param} key={key} onChange={this.props.onChange} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} name={key} value={values[key]} />
            })}
        </div>;
    }
}

function controlFor(param: PD.Any): ParamControl | undefined {
    switch (param.type) {
        case 'value': return void 0;
        case 'boolean': return BoolControl;
        case 'number': return NumberControl;
        case 'converted': return ConvertedControl;
        case 'multi-select': return MultiSelectControl;
        case 'color': return ColorControl;
        case 'select': return SelectControl;
        case 'text': return TextControl;
        case 'interval': return IntervalControl;
        case 'group': return GroupControl;
        case 'mapped': return MappedControl;
        case 'line-graph': return void 0;
    }
    throw new Error('not supported');
}

// type ParamWrapperProps = { name: string, value: any, param: PD.Base<any>, onChange: ParamOnChange, control: ValueControl, onEnter?: () => void, isEnabled?: boolean }

export type ParamOnChange = (params: { param: PD.Base<any>, name: string, value: any }) => void
export interface ParamProps<P extends PD.Base<any> = PD.Base<any>> { name: string, value: P['defaultValue'], param: P, isDisabled?: boolean, onChange: ParamOnChange, onEnter?: () => void }
export type ParamControl = React.ComponentClass<ParamProps<any>>

export abstract class SimpleParam<P extends PD.Any> extends React.PureComponent<ParamProps<P>> {
    protected update(value: any) {
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    abstract renderControl(): JSX.Element;

    render() {
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <div style={{ padding: '0 3px', borderBottom: '1px solid #ccc' }}>
            <div style={{ lineHeight: '20px', float: 'left' }} title={this.props.param.description}>{label}</div>
            <div style={{ float: 'left', marginLeft: '5px' }}>
                {this.renderControl()}
            </div>
            <div style={{ clear: 'both' }} />
        </div>;
    }
}

export class BoolControl extends SimpleParam<PD.Boolean> {
    onClick = () => { this.update(!this.props.value); }
    renderControl() {
        return <button onClick={this.onClick} disabled={this.props.isDisabled}>{this.props.value ? '✓ On' : '✗ Off'}</button>;
    }
}

export class NumberControl extends SimpleParam<PD.Numeric> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => { this.update(+e.target.value); }
    renderControl() {
        return <span>
            <input type='range' value={'' + this.props.value} min={this.props.param.min} max={this.props.param.max} step={this.props.param.step} onChange={this.onChange} disabled={this.props.isDisabled} />
            <br />{this.props.value}
        </span>
    }
}

export class TextControl extends SimpleParam<PD.Text> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const value = e.target.value;
        if (value !== this.props.value) {
            this.update(value);
        }
    }

    onKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if (!this.props.onEnter) return;
        if ((e.keyCode === 13 || e.charCode === 13)) {
            this.props.onEnter();
        }
    }

    renderControl() {
        return <input type='text'
            value={this.props.value || ''}
            onChange={this.onChange}
            onKeyPress={this.props.onEnter ? this.onKeyPress : void 0}
            disabled={this.props.isDisabled}
        />;
    }
}

export class SelectControl extends SimpleParam<PD.Select<any>> {
    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => { this.update(e.target.value); }
    renderControl() {
        return <select value={this.props.value || ''} onChange={this.onChange} disabled={this.props.isDisabled}>
            {this.props.param.options.map(([value, label]) => <option key={value} value={value}>{label}</option>)}
        </select>;
    }
}

export class IntervalControl extends SimpleParam<PD.Interval> {
    // onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    //     this.setState({ value: e.target.value });
    //     this.props.onChange(e.target.value);
    // }

    renderControl() {
        return <span>interval TODO</span>;
    }
}

export class ColorControl extends SimpleParam<PD.Color> {
    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        this.update(Color(parseInt(e.target.value)));
    }

    renderControl() {
        return <select value={this.props.value} onChange={this.onChange}>
            {Object.keys(ColorNames).map(name => {
                return <option key={name} value={(ColorNames as { [k: string]: Color})[name]}>{name}</option>
            })}
        </select>;
    }
}

export class MultiSelectControl extends React.PureComponent<ParamProps<PD.MultiSelect<any>>> {
    change(value: PD.MultiSelect<any>['defaultValue'] ) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    toggle(key: string) {
        return () => {
            if (this.props.value.indexOf(key) < 0) this.change(this.props.value.concat(key));
            else this.change(this.props.value.filter(v => v !== key))
        }
    }

    render() {
        const current = this.props.value;
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <div>
            <div>{label} <small>{`${current.length} of ${this.props.param.options.length}`}</small></div>
            <div style={{ paddingLeft: '7px' }}>
                {this.props.param.options.map(([value, label]) =>
                    <button key={value} onClick={this.toggle(value)} disabled={this.props.isDisabled}>
                        {current.indexOf(value) >= 0 ? `✓ ${label}` : `✗ ${label}`}
                    </button>)}
            </div>
        </div>;
    }
}

export class GroupControl extends React.PureComponent<ParamProps<PD.Group<any>>> {
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
        const label = this.props.param.label || camelCaseToWords(this.props.name);

        // TODO toggle panel
        return <div>
            <div>{label}</div>
            <ParameterControls params={params} onChange={this.onChangeParam} values={value.params} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />
        </div>
    }
}

export class MappedControl extends React.PureComponent<ParamProps<PD.Mapped<any>>> {
    change(value: PD.Mapped<any>['defaultValue'] ) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeName: ParamOnChange = e => {
        // TODO: Cache values when changing types?
        this.change({ name: e.value, params: this.props.param.map(e.value).defaultValue });
    }

    onChangeParam: ParamOnChange = e => {
        const value: PD.Mapped<any>['defaultValue'] = this.props.value;
        this.change({ name: value.name, params: e.value });
    }

    render() {
        const value: PD.Mapped<any>['defaultValue'] = this.props.value;
        const param = this.props.param.map(value.name);
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const Mapped = controlFor(param);

        const select = <SelectControl param={this.props.param.select}
            isDisabled={this.props.isDisabled} onChange={this.onChangeName} onEnter={this.props.onEnter}
            name={label} value={value.name} />

        if (!Mapped) {
            return select;
        }

        return <div>
            {select}
            <div style={{ borderLeft: '5px solid #777', paddingLeft: '5px' }}>
                <Mapped param={param} value={value} name='' onChange={this.onChangeParam} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />
            </div>
        </div>
    }
}

export class ConvertedControl extends React.PureComponent<ParamProps<PD.Converted<any, any>>> {
    onChange: ParamOnChange = e => {
        this.props.onChange({
            name: this.props.name,
            param: this.props.param,
            value: this.props.param.toValue(e.value)
        });
    }

    render() {
        const value = this.props.param.fromValue(this.props.value);
        const Converted = controlFor(this.props.param.converted);

        if (!Converted) return null;
        return <Converted param={this.props.param.converted} value={value} name={this.props.name} onChange={this.onChange} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />
    }
}
