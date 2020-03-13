/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Color } from '../../mol-util/color';
import { PurePluginUIComponent } from '../base';
import { IconName, Icon } from './icons';
import { ShapeName, Shape } from './shapes';

export class ControlGroup extends React.Component<{
    header: string,
    initialExpanded?: boolean,
    hideExpander?: boolean,
    hideOffset?: boolean,
    topRightIcon?: IconName,
    onHeaderClick?: () => void
}, { isExpanded: boolean }> {
    state = { isExpanded: !!this.props.initialExpanded }

    headerClicked = () => {
        if (this.props.onHeaderClick) {
            this.props.onHeaderClick();
        } else {
            this.setState({ isExpanded: !this.state.isExpanded });
        }
    }

    render() {
        // TODO: customize header style (bg color, togle button etc)
        return <div className='msp-control-group-wrapper' style={{ position: 'relative' }}>
            <div className='msp-control-group-header'>
                <button className='msp-btn msp-btn-block' onClick={this.headerClicked}>
                    {!this.props.hideExpander && <Icon name={this.state.isExpanded ? 'collapse' : 'expand'} />}
                    {this.props.topRightIcon && <Icon name={this.props.topRightIcon} style={{ position: 'absolute', right: '2px', top: 0 }} />}
                    {this.props.header}
                </button>
            </div>
            {this.state.isExpanded && <div className={this.props.hideOffset ? '' : 'msp-control-offset'} style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                {this.props.children}
            </div>
            }
        </div>
    }
}

export interface TextInputProps<T> {
    className?: string,
    style?: React.CSSProperties,
    value: T,
    fromValue?(v: T): string,
    toValue?(s: string): T,
    // TODO: add error/help messages here?
    isValid?(s: string): boolean,
    onChange(value: T): void,
    onEnter?(): void,
    onBlur?(): void,
    delayMs?: number,
    blurOnEnter?: boolean,
    blurOnEscape?: boolean,
    isDisabled?: boolean,
    placeholder?: string
}

interface TextInputState {
    originalValue: string,
    value: string
}

function _id(x: any) { return x; }

export class TextInput<T = string> extends PurePluginUIComponent<TextInputProps<T>, TextInputState> {
    private input = React.createRef<HTMLInputElement>();
    private delayHandle: any = void 0;
    private pendingValue: T | undefined = void 0;

    state = { originalValue: '', value: '' }

    onBlur = () => {
        this.setState({ value: '' + this.state.originalValue });
        if (this.props.onBlur) this.props.onBlur();
    }

    get isPending() { return typeof this.delayHandle !== 'undefined'; }

    clearTimeout() {
        if (this.isPending) {
            clearTimeout(this.delayHandle);
            this.delayHandle = void 0;
        }
    }

    raiseOnChange = () => {
        this.props.onChange(this.pendingValue!);
        this.pendingValue = void 0;
    }

    triggerChanged(formatted: string, converted: T) {
        this.clearTimeout();

        if (formatted === this.state.originalValue) return;

        if (this.props.delayMs) {
            this.pendingValue = converted;
            this.delayHandle = setTimeout(this.raiseOnChange, this.props.delayMs);
        } else {
            this.props.onChange(converted);
        }
    }

    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const value = e.target.value;

        if (this.props.isValid && !this.props.isValid(value)) {
            this.clearTimeout();
            this.setState({ value });
            return;
        }

        const converted = (this.props.toValue || _id)(value);
        const formatted = (this.props.fromValue || _id)(converted);
        this.setState({ value: formatted }, () => this.triggerChanged(formatted, converted));
    }

    onKeyUp  = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if (e.charCode === 27 || e.keyCode === 27 /* esc */) {
            if (this.props.blurOnEscape && this.input.current) {
                this.input.current.blur();
            }
        }
    }

    onKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if (e.keyCode === 13 || e.charCode === 13 /* enter */) {
            if (this.isPending) {
                this.clearTimeout();
                this.raiseOnChange();
            }
            if (this.props.blurOnEnter && this.input.current) {
                this.input.current.blur();
            }
            if (this.props.onEnter) this.props.onEnter();
        }
        e.stopPropagation();
    }

    static getDerivedStateFromProps(props: TextInputProps<any>, state: TextInputState) {
        const value = props.fromValue ? props.fromValue(props.value) : props.value;
        if (value === state.originalValue) return null;
        return { originalValue: value, value };
    }

    render() {
        return <input type='text'
            className={this.props.className}
            style={this.props.style}
            ref={this.input}
            onBlur={this.onBlur}
            value={this.state.value}
            placeholder={this.props.placeholder}
            onChange={this.onChange}
            onKeyPress={this.props.onEnter || this.props.blurOnEnter || this.props.blurOnEscape ? this.onKeyPress : void 0}
            onKeyDown={this.props.blurOnEscape ? this.onKeyUp : void 0}
            disabled={!!this.props.isDisabled}
        />;
    }
}

// TODO: replace this with parametrized TextInput
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
        e.stopPropagation();
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
                        <Icon name={this.state.isExpanded ? 'minus' : 'plus'} style={{ display: 'inline-block' }} />
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

export function IconButton(props: {
    icon: IconName,
    small?: boolean,
    onClick: (e: React.MouseEvent<HTMLButtonElement>) => void,
    title?: string,
    toggleState?: boolean,
    disabled?: boolean,
    customClass?: string,
    style?: React.CSSProperties,
    'data-id'?: string,
    extraContent?: JSX.Element
}) {
    let className = `msp-btn-link msp-btn-icon${props.small ? '-small' : ''}${props.customClass ? ' ' + props.customClass : ''}`;
    if (typeof props.toggleState !== 'undefined') {
        className += ` msp-btn-link-toggle-${props.toggleState ? 'on' : 'off'}`
    }
    const iconStyle = props.small ? { fontSize: '80%' } : void 0;
    return <button className={className} onClick={props.onClick} title={props.title} disabled={props.disabled} data-id={props['data-id']} style={props.style}>
        <Icon name={props.icon} style={iconStyle} />
        {props.extraContent}
    </button>;
}

export class ButtonSelect extends React.PureComponent<{ label: string, onChange: (value: string) => void, disabled?: boolean }> {
    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        e.preventDefault()
        this.props.onChange(e.target.value)
        e.target.value = '_'
    }

    render() {
        return <select value='_' onChange={this.onChange} disabled={this.props.disabled}>
            <option key='_' value='_'>{this.props.label}</option>
            {this.props.children}
        </select>
    }
}

export function Options(options: [string, string][]) {
    return options.map(([value, label]) => <option key={value} value={value}>{label}</option>)
}

export function SectionHeader(props: { icon?: IconName, title: string | JSX.Element, desc?: string}) {
    return <div className='msp-section-header'>
        {props.icon && <Icon name={props.icon} />}
        {props.title} <small>{props.desc}</small>
    </div>
}

export type ToggleButtonProps = {
    style?: React.CSSProperties,
    className?: string,
    disabled?: boolean,
    label?: string | JSX.Element,
    title?: string,
    icon?: IconName,
    shape?: ShapeName,
    isSelected?: boolean,
    toggle: () => void
}

export class ToggleButton extends React.PureComponent<ToggleButtonProps> {
    onClick = (e: React.MouseEvent<HTMLButtonElement>) => {
        e.currentTarget.blur();
        this.props.toggle();
    }

    render() {
        const props = this.props;
        const label = props.label;
        return <button onClick={this.onClick} title={this.props.title}
            disabled={props.disabled} style={props.style} className={props.className}>
            <Icon name={this.props.icon} />
            <Shape name={this.props.shape} />
            {label && this.props.isSelected ? <b>{label}</b> : label}
        </button>;
    }
}

export class ExpandGroup extends React.PureComponent<{ header: string, initiallyExpanded?: boolean, noOffset?: boolean, marginTop?: 0 | string }, { isExpanded: boolean }> {
    state = { isExpanded: !!this.props.initiallyExpanded };

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    render() {
        return <>
            <div className='msp-control-group-header' style={{ marginTop: this.props.marginTop !== void 0 ? this.props.marginTop : '1px' }}>
                <button className='msp-btn msp-btn-block' onClick={this.toggleExpanded}>
                    <Icon name={this.state.isExpanded ? 'collapse' : 'expand'} />
                    {this.props.header}
                </button>
            </div>
            {this.state.isExpanded &&
                (this.props.noOffset
                    ? this.props.children
                    : <div className='msp-control-offset'>
                        {this.props.children}
                    </div>)}
        </>;
    }
}