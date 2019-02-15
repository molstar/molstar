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