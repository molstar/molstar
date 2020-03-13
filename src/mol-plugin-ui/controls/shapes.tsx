/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';

export type ShapeName =
    | '' | 'set' | 'intersect' | 'union' | 'subtract'

export function Shape(props: {
    name?: ShapeName,
    style?: React.CSSProperties,
    title?: string
}) {
    const shape = getShape(props.name)
    if (!shape) return null
    return <span style={props.style} title={props.title}>
        {shape}
    </span>;
}

function getShape(name?: ShapeName) {
    switch (name) {
        case 'union': return <Union />
        case 'subtract': return <Subtract />
        case 'intersect': return <Intersect />
        case 'set': return <Set />
        default: return null
    }
}

const circleLeft = <circle r="6px" id="circle-left" cy="16px" cx="12px" strokeWidth="1"/>
const circleRight = <circle r="6px" id="circle-right" cy="16px" cx="20px" strokeWidth="1"/>

function Union() {
    return <svg width="32px" height="32px">
        <defs>
            {circleLeft}
            {circleRight}
        </defs>
        <g>
            <use href="#circle-left" className="msp-shape-filled"/>
            <use href="#circle-right" className="msp-shape-filled"/>
        </g>
    </svg>;
}

function Subtract() {
    return <svg width="32px" height="32px">
        <defs>
            {circleLeft}
            {circleRight}
            <mask id="mask-left">
                <use href="#circle-left" fill="white" stroke="white"/>
                <use href="#circle-right" fill="black" strokeWidth="0" stroke="white"/>
            </mask>
        </defs>
        <g>
            <use href="#circle-left" className="msp-shape-filled" mask="url(#mask-left)"/>
            <use href="#circle-right" className="msp-shape-empty"/>
        </g>
    </svg>;
}

function Intersect() {
    return <svg width="32px" height="32px">
        <defs>
            {circleLeft}
            {circleRight}
            <clipPath id="clip-left">
                <use href="#circle-right"/>
            </clipPath>
        </defs>
        <g>
            <use href="#circle-left" className="msp-shape-filled" clipPath="url(#clip-left)"/>
            <use href="#circle-left" className="msp-shape-empty"/>
            <use href="#circle-right" className="msp-shape-empty"/>
        </g>
    </svg>;
}

function Set() {
    return <svg width="32px" height="32px">
        <defs>
            {circleLeft}
            {circleRight}
        </defs>
        <g>
            <use href="#circle-left" className="msp-shape-filled"/>
            <use href="#circle-right" className="msp-shape-empty"/>
        </g>
    </svg>;
}