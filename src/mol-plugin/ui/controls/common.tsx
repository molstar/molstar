/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';

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