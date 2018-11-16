/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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