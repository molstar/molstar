/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { camelCaseToWords } from 'mol-util/string';
import { ColorNames, ColorNamesValueMap } from 'mol-util/color/tables';
import { Color } from 'mol-util/color';
import { Vec2 } from 'mol-math/linear-algebra';
import LineGraphComponent from './line-graph/line-graph-component';

import { Slider, Slider2 } from './slider';


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
                if (param.isHidden) return null;
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
        case 'number': return typeof param.min !== 'undefined' && typeof param.max !== 'undefined'
            ? NumberRangeControl : NumberInputControl;
        case 'converted': return ConvertedControl;
        case 'multi-select': return MultiSelectControl;
        case 'color': return ColorControl;
        case 'vec3': return Vec3Control;
        case 'file': return FileControl;
        case 'select': return SelectControl;
        case 'text': return TextControl;
        case 'interval': return typeof param.min !== 'undefined' && typeof param.max !== 'undefined'
        ? BoundedIntervalControl : IntervalControl;
        case 'group': return GroupControl;
        case 'mapped': return MappedControl;
        case 'line-graph': return LineGraphControl;
    }
    console.warn(`${(param as any).type} has no associated UI component.`);
    return void 0;
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
        return <div className='msp-control-row'>
            <span title={this.props.param.description}>{label}</span>
            <div>
                {this.renderControl()}
            </div>
        </div>;
    }
}

export class BoolControl extends SimpleParam<PD.Boolean> {
    onClick = (e: React.MouseEvent<HTMLButtonElement>) => { this.update(!this.props.value); e.currentTarget.blur(); }
    renderControl() {
        return <button onClick={this.onClick} disabled={this.props.isDisabled}>
            <span className={`msp-icon msp-icon-${this.props.value ? 'ok' : 'off'}`} />
            {this.props.value ? 'On' : 'Off'}
        </button>;
    }
}

export class LineGraphControl extends React.PureComponent<ParamProps<PD.LineGraph>, { isExpanded: boolean, isOverPoint: boolean, message: string }> {
    state = { 
        isExpanded: false,
        isOverPoint: false,
        message: `${this.props.param.defaultValue.length} points`,
    }

    onHover = (point?: Vec2) => {
        this.setState({isOverPoint: !this.state.isOverPoint});
        if(point){
            this.setState({message: `(${point[0].toFixed(2)}, ${point[1].toFixed(2)})`});
            return;
        }
        this.setState({message: `${this.props.value.length} points`});
    }

    onDrag = (point: Vec2) => {
        this.setState({message: `(${point[0].toFixed(2)}, ${point[1].toFixed(2)})`});
    }

    onChange = (value: PD.LineGraph['defaultValue'] ) => {
        this.props.onChange({ name: this.props.name, param: this.props.param, value: value});
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    render() {
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <>
            <div className='msp-control-row'>
                <span>{label}</span>
                <div>
                    <button onClick={this.toggleExpanded}>
                        {`${this.state.message}`}
                    </button>
                </div>
            </div>
            <div className='msp-control-offset' style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                <LineGraphComponent
                    data={this.props.param.defaultValue} 
                    onChange={this.onChange} 
                    onHover={this.onHover}
                    onDrag={this.onDrag}/>
            </div>
        </>;
    }
}

export class NumberInputControl extends SimpleParam<PD.Numeric> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => { this.update(+e.target.value); }
    renderControl() {
        return <span>
            number input TODO
        </span>
    }
}

export class NumberRangeControl extends SimpleParam<PD.Numeric> {
    onChange = (v: number) => { this.update(v); }
    renderControl() {
        return <Slider value={this.props.value} min={this.props.param.min!} max={this.props.param.max!}
            step={this.props.param.step} onChange={this.onChange} disabled={this.props.isDisabled} />
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
        const placeholder = this.props.param.label || camelCaseToWords(this.props.name);
        return <input type='text'
            value={this.props.value || ''}
            placeholder={placeholder}
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
    onChange = (v: [number, number]) => { this.update(v); }
    renderControl() {
        return <span>interval TODO</span>;
    }
}

export class BoundedIntervalControl extends SimpleParam<PD.Interval> {
    onChange = (v: [number, number]) => { this.update(v); }
    renderControl() {
        return <Slider2 value={this.props.value} min={this.props.param.min!} max={this.props.param.max!}
            step={this.props.param.step} onChange={this.onChange} disabled={this.props.isDisabled} />;
    }
}

let _colors: React.ReactFragment | undefined = void 0;
function ColorOptions() {
    if (_colors) return _colors;
    _colors = <>{Object.keys(ColorNames).map(name =>
        <option key={name} value={(ColorNames as { [k: string]: Color})[name]} style={{ background: `${Color.toStyle((ColorNames as { [k: string]: Color})[name])}` }} >
            {name}
        </option>
    )}</>;
    return _colors;
}

function ColorValueOption(color: Color) {
    return !ColorNamesValueMap.has(color) ? <option key={Color.toHexString(color)} value={color} style={{ background: `${Color.toStyle(color)}` }} >
        {Color.toHexString(color)}
    </option> : null
}


export class ColorControl extends SimpleParam<PD.Color> {
    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        this.update(Color(parseInt(e.target.value)));
    }

    renderControl() {
        return <select value={this.props.value} onChange={this.onChange} style={{ borderLeft: `16px solid ${Color.toStyle(this.props.value)}` }}>
            {ColorValueOption(this.props.value)}
            {ColorOptions()}
        </select>;
    }
}

export class Vec3Control extends SimpleParam<PD.Vec3> {
    // onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    //     this.setState({ value: e.target.value });
    //     this.props.onChange(e.target.value);
    // }

    renderControl() {
        return <span>vec3 TODO</span>;
    }
}

export class FileControl extends React.PureComponent<ParamProps<PD.FileParam>> {
    change(value: File) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeFile = (e: React.ChangeEvent<HTMLInputElement>) => {
        this.change(e.target.files![0]);
    }

    render() {
        const value = this.props.value;

        // return <input disabled={this.props.isDisabled} value={void 0} type='file' multiple={false} />
        return <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file' style={{ marginTop: '1px' }}>
            {value ? value.name : 'Select a file...'} <input disabled={this.props.isDisabled} onChange={this.onChangeFile} type='file' multiple={false} accept={this.props.param.accept} />
        </div>
    }
}

export class MultiSelectControl extends React.PureComponent<ParamProps<PD.MultiSelect<any>>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    change(value: PD.MultiSelect<any>['defaultValue'] ) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    toggle(key: string) {
        return (e: React.MouseEvent<HTMLButtonElement>) => {
            if (this.props.value.indexOf(key) < 0) this.change(this.props.value.concat(key));
            else this.change(this.props.value.filter(v => v !== key));
            e.currentTarget.blur();
        }
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    render() {
        const current = this.props.value;
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <>
            <div className='msp-control-row'>
                <span>{label}</span>
                <div>
                    <button onClick={this.toggleExpanded}>
                        {`${current.length} of ${this.props.param.options.length}`}
                    </button>
                </div>
            </div>
            <div className='msp-control-offset' style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                {this.props.param.options.map(([value, label]) => {
                    const sel = current.indexOf(value) >= 0;
                    return <div key={value} className='msp-row'>
                        <button onClick={this.toggle(value)} disabled={this.props.isDisabled}>
                            <span style={{ float: sel ? 'left' : 'right' }}>{sel ? `✓ ${label}` : `${label} ✗`}</span>
                        </button>
                </div> })}
            </div>
        </>;
    }
}

export class GroupControl extends React.PureComponent<ParamProps<PD.Group<any>>, { isExpanded: boolean }> {
    state = { isExpanded: !!this.props.param.isExpanded }

    change(value: any ) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeParam: ParamOnChange = e => {
        this.change({ ...this.props.value, [e.name]: e.value });
    }

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    render() {
        const params = this.props.param.params;
        const label = this.props.param.label || camelCaseToWords(this.props.name);

        const controls = <ParameterControls params={params} onChange={this.onChangeParam} values={this.props.value} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />;

        if (this.props.param.isFlat) {
            return controls;
        }

        return <div className='msp-control-group-wrapper'>
            <div className='msp-control-group-header'>
                <button className='msp-btn msp-btn-block' onClick={this.toggleExpanded}>
                    <span className={`msp-icon msp-icon-${this.state.isExpanded ? 'collapse' : 'expand'}`} />
                    {label}
                </button>
            </div>
            {this.state.isExpanded && <div className='msp-control-offset' style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                {controls}
            </div>
            }
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
        this.change({ name: this.props.value.name, params: e.value });
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
            <Mapped param={param} value={value.params} name={`${label} Properties`} onChange={this.onChangeParam} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />
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
