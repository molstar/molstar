/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';

export function Icon(props: {
    svg?: React.FC,
    // name?: IconName | undefined,
    style?: React.CSSProperties,
    title?: string
}) {
    if (!props.svg) return null;
    return <span className='msp-icon msp-material-icon' title={props.title} style={props.style}><props.svg /></span>;
}

//

export function Union() { return _union; }
export function Subtract() { return _subtract; }
export function Intersect() { return _intersect; }
export function SetSvg() { return _set; }

const circleLeft = <circle r="6px" id="circle-left" cy="12px" cx="8px" strokeWidth="0.5" />;
const circleRight = <circle r="6px" id="circle-right" cy="12px" cx="16px" strokeWidth="0.5" />;

const _union = <svg width="24px" height="24px" viewBox="0 0 24 24">
    <defs>
        {circleLeft}
        {circleRight}
    </defs>
    <g>
        <use href="#circle-left" className="msp-shape-filled" />
        <use href="#circle-right" className="msp-shape-filled" />
        <use href="#circle-left" className="msp-shape-empty" />
    </g>
</svg>;

const _subtract = <svg width="24px" height="24px" viewBox="0 0 24 24">
    <defs>
        {circleLeft}
        {circleRight}
        <mask id="mask-left">
            <use href="#circle-left" fill="white" stroke="white" />
            <use href="#circle-right" fill="black" strokeWidth="0" stroke="white" />
        </mask>
    </defs>
    <g>
        <use href="#circle-left" className="msp-shape-filled" mask="url(#mask-left)" />
        <use href="#circle-right" className="msp-shape-empty" />
    </g>
</svg>;

const _intersect = <svg width="24px" height="24px" viewBox="0 0 24 24">
    <defs>
        {circleLeft}
        {circleRight}
        <clipPath id="clip-left">
            <use href="#circle-right" />
        </clipPath>
    </defs>
    <g>
        <use href="#circle-left" className="msp-shape-filled" clipPath="url(#clip-left)" />
        <use href="#circle-left" className="msp-shape-empty" />
        <use href="#circle-right" className="msp-shape-empty" />
    </g>
</svg>;

const _set = <svg width="24px" height="24px" viewBox="0 0 24 24">
    <defs>
        {circleLeft}
        {circleRight}
    </defs>
    <g>
        <use href="#circle-left" className="msp-shape-empty" />
        <use href="#circle-right" className="msp-shape-filled" />
    </g>
</svg>;