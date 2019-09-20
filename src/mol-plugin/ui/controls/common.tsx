/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Color } from '../../../mol-util/color';

export class ControlGroup extends React.Component<{ header: string, initialExpanded?: boolean }, { isExpanded: boolean }> {
    state = { isExpanded: !!this.props.initialExpanded }

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    render() {
        return <div className='msp-control-group-wrapper'>
            <div className='msp-control-group-header'>
                <button className='msp-btn msp-btn-block' onClick={this.toggleExpanded}>
                    <span className={`msp-icon msp-icon-${this.state.isExpanded ? 'collapse' : 'expand'}`} />
                    {this.props.header}
                </button>
            </div>
            {this.state.isExpanded && <div className='msp-control-offset' style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                {this.props.children}
            </div>
            }
        </div>
    }
}

export class NumericInput extends React.PureComponent<{
    value: number,
    onChange: (v: number) => void,
    onEnter?: () => void,
    onBlur?: () => void,
    blurOnEnter?: boolean,
    isDisabled?: boolean,
    placeholder?: string
}, { value: string }> {
    state = { value: '0' };
    input = React.createRef<HTMLInputElement>();

    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const value = +e.target.value;
        this.setState({ value: e.target.value }, () => {
            if (!Number.isNaN(value) && value !== this.props.value) {
                this.props.onChange(value);
            }
        });
    }

    onKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if ((e.keyCode === 13 || e.charCode === 13)) {
            if (this.props.blurOnEnter && this.input.current) {
                this.input.current.blur();
            }
            if (this.props.onEnter) this.props.onEnter();
        }
    }

    onBlur = () => {
        this.setState({ value: '' + this.props.value });
        if (this.props.onBlur) this.props.onBlur();
    }

    static getDerivedStateFromProps(props: { value: number }, state: { value: string }) {
        const value = +state.value;
        if (Number.isNaN(value) || value === props.value) return null;
        return { value: '' + props.value };
    }

    render() {
        return <input type='text'
            ref={this.input}
            onBlur={this.onBlur}
            value={this.state.value}
            placeholder={this.props.placeholder}
            onChange={this.onChange}
            onKeyPress={this.props.onEnter || this.props.blurOnEnter ? this.onKeyPress : void 0}
            disabled={!!this.props.isDisabled}
        />
    }
}

export function Icon(props: {
    name: string
}) {
    return <span className={`msp-icon msp-icon-${props.name}`} />;
}

export function IconButton(props: {
    icon: string,
    isSmall?: boolean,
    onClick: (e: React.MouseEvent<HTMLButtonElement>) => void,
    title?: string,
    toggleState?: boolean,
    disabled?: boolean,
    customClass?: string,
    'data-id'?: string
}) {
    let className = `msp-btn-link msp-btn-icon${props.isSmall ? '-small' : ''}${props.customClass ? ' ' + props.customClass : ''}`;
    if (typeof props.toggleState !== 'undefined') className += ` msp-btn-link-toggle-${props.toggleState ? 'on' : 'off'}`
    return <button className={className} onClick={props.onClick} title={props.title} disabled={props.disabled} data-id={props['data-id']}>
        <span className={`msp-icon msp-icon-${props.icon}`}/>
    </button>;
}

export class ExpandableGroup extends React.Component<{
    label: string,
    colorStripe?: Color,
    pivot: JSX.Element,
    controls: JSX.Element
}, { isExpanded: boolean }> {
    state = { isExpanded: false };

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    render() {
        const { label, pivot, controls } = this.props;
        // TODO: fix the inline CSS
        return <>
            <div className='msp-control-row'>
                <span>
                    {label}
                    <button className='msp-btn-link msp-btn-icon msp-control-group-expander' onClick={this.toggleExpanded} title={`${this.state.isExpanded ? 'Less' : 'More'} options`}
                        style={{ background: 'transparent', textAlign: 'left', padding: '0' }}>
                        <span className={`msp-icon msp-icon-${this.state.isExpanded ? 'minus' : 'plus'}`} style={{ display: 'inline-block' }} />
                    </button>
                </span>
                <div>{pivot}</div>
                {this.props.colorStripe && <div className='msp-expandable-group-color-stripe' style={{ backgroundColor: Color.toStyle(this.props.colorStripe) }} /> }
            </div>
            {this.state.isExpanded && <div className='msp-control-offset'>
                {controls}
            </div>}
        </>;
    }
}

export class ButtonSelect extends React.PureComponent<{ label: string, onChange: (value: string) => void }> {

    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        e.preventDefault()
        this.props.onChange(e.target.value)
        e.target.value = '_'
    }

    render() {
        return <select value='_' onChange={this.onChange}>
            <option key='_' value='_'>{this.props.label}</option>
            {this.props.children}
        </select>
    }
}

export function Options(options: [string, string][]) {
    return options.map(([value, label]) => <option key={value} value={value}>{label}</option>)
}

// export const ToggleButton = (props: {
//     onChange: (v: boolean) => void,
//     value: boolean,
//     label: string,
//     title?: string
// }) => <div className='lm-control-row lm-toggle-button' title={props.title}>
//         <span>{props.label}</span>
//         <div>
//             <button onClick={e => { props.onChange.call(null, !props.value); (e.target as HTMLElement).blur(); }}>
//                     <span className={ `lm-icon lm-icon-${props.value ? 'ok' : 'off'}` }></span> {props.value ? 'On' : 'Off'}
//             </button>
//         </div>
//     </div>