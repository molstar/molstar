/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Subject } from 'rxjs';

export function createChangeSubject(): ParamChanges {
    return new Subject<{ param: PD.Base<any>, name: string, value: any }>();
}

export interface ParameterControlsProps<P extends PD.Params = PD.Params> {
    params: P,
    values: any,
    changes: ParamChanges,
    isEnabled?: boolean,
    onEnter?: () => void
}

export class ParameterControls<P extends PD.Params> extends React.PureComponent<ParameterControlsProps<P>, {}> {
    render() {
        const common = {
            changes: this.props.changes,
            isEnabled: this.props.isEnabled,
            onEnter: this.props.onEnter,
        }
        const params = this.props.params;
        const values = this.props.values;
        return <div style={{ width: '100%' }}>
            {Object.keys(params).map(key => <ParamWrapper control={controlFor(params[key])} param={params[key]} key={key} {...common} name={key} value={values[key]} />)}
        </div>;
    }
}

function controlFor(param: PD.Any): ValueControl {
    switch (param.type) {
        case 'boolean': return BoolControl;
        case 'number': return NumberControl;
        case 'range': return NumberControl;
        case 'multi-select': throw new Error('nyi');
        case 'color': throw new Error('nyi');
        case 'select': return SelectControl;
        case 'text': return TextControl;
    }
    throw new Error('not supporter');
}
type ParamWrapperProps = { name: string, value: any, param: PD.Base<any>, changes: ParamChanges, control: ValueControl, onEnter?: () => void, isEnabled?: boolean }
export type ParamChanges = Subject<{ param: PD.Base<any>, name: string, value: any }>
type ValueControlProps<P extends PD.Base<any> = PD.Base<any>> = { value: any, param: P, isEnabled?: boolean, onChange: (v: any) => void, onEnter?: () => void }
type ValueControl = React.ComponentClass<ValueControlProps<any>>

export class ParamWrapper extends React.PureComponent<ParamWrapperProps> {
    onChange = (value: any) => {
        this.props.changes.next({ param: this.props.param, name: this.props.name, value });
    }

    render() {
        return <div>
            <span title={this.props.param.description}>{this.props.param.label}</span>
            <div>
                <this.props.control value={this.props.value} param={this.props.param} onChange={this.onChange} onEnter={this.props.onEnter} isEnabled={this.props.isEnabled} />
            </div>
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

export class NumberControl extends React.PureComponent<ValueControlProps<PD.Numeric>> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        this.props.onChange(+e.target.value);
    }

    render() {
        return <input type='range'
            value={this.props.value}
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
        return <select value={this.props.value} onChange={this.onChange}>
            {this.props.param.options.map(([value, label]) => <option key={label} value={value}>{label}</option>)}
        </select>;
    }
}