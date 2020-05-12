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
    title?: string,
    /** Adds right margin to the icon */
    inline?: boolean
}) {
    if (!props.svg) return null;
    return <span className={`msp-icon msp-material-icon${props.inline ? ' msp-icon-inline' : ''}`} title={props.title} style={props.style}><props.svg /></span>;
}

//

const circleLeft = <circle r="6px" id="circle-left" cy="12px" cx="8px" strokeWidth="1" />;
const circleRight = <circle r="6px" id="circle-right" cy="12px" cx="16px" strokeWidth="1" />;

export function UnionSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
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
}

export function SubtractSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
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
}

export function IntersectSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
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
}

export function SetSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <defs>
            {circleLeft}
            {circleRight}
        </defs>
        <g>
            <use href="#circle-left" className="msp-shape-empty" />
            <use href="#circle-right" className="msp-shape-filled" />
        </g>
    </svg>;
}

export function MoleculeSvg() {
    return <svg width="17px" height="17px" viewBox="0 0 299.463 299.463">
        <g>
            <path d="M256.851,173.832v-48.201c22.916-4.918,34.151-30.668,22.556-50.771c-11.547-20.004-39.486-23.251-55.242-5.844 l-41.746-24.106C189.618,22.603,172.861,0,149.734,0c-23.132,0-39.881,22.609-32.685,44.911L75.305,69.016 C59.522,51.586,31.597,54.88,20.061,74.863c-11.63,20.163-0.298,45.862,22.557,50.769v48.2 c-22.821,4.898-34.195,30.591-22.556,50.771c11.529,19.972,39.454,23.285,55.242,5.845l41.746,24.106 c-7.199,22.308,9.559,44.911,32.685,44.911c23.132,0,39.88-22.609,32.685-44.911l41.745-24.106 c15.817,17.469,43.73,14.099,55.242-5.844c0,0,0-0.001,0.001-0.002c4.587-7.953,5.805-17.213,3.431-26.076 C279.392,185.657,269.129,176.461,256.851,173.832z M249.62,72.088c20.568,0,27.428,27.191,10.008,37.239 c-0.003,0.002-0.006,0.003-0.009,0.005c-10.04,5.81-22.85,1.762-27.877-8.475C225.206,87.548,234.938,72.088,249.62,72.088z M149.734,14.4c11.005,0,19.958,8.954,19.958,19.959c0,11.127-9.077,19.958-19.958,19.958c-10.95,0-19.958-8.9-19.958-19.958 C129.776,23.354,138.729,14.4,149.734,14.4z M39.84,109.328c-17.451-10.067-10.534-37.24,10.01-37.24 c15.311,0,24.922,16.653,17.251,29.942C61.681,111.397,49.517,114.925,39.84,109.328z M59.802,224.702 c-9.535,5.503-21.768,2.229-27.268-7.298c-7.639-13.242,1.887-29.945,17.236-29.945c0.013,0,0.027,0,0.04,0 C70.07,187.48,77.49,214.469,59.802,224.702z M149.734,285.062c-11.005,0-19.958-8.954-19.958-19.958 c0-11.127,9.077-19.958,19.958-19.958c10.954,0,19.958,8.903,19.958,19.958C169.693,276.109,160.74,285.062,149.734,285.062z M216.953,217.982l-41.727,24.095c-13.778-15.22-37.459-14.94-50.983,0l-41.728-24.096c6.196-19.289-5.541-39.835-25.498-44.149 V125.63c19.752-4.268,31.762-24.65,25.498-44.149l41.727-24.095c13.629,15.055,37.32,15.093,50.983,0l41.728,24.096 c-6.196,19.29,5.534,39.835,25.498,44.149v48.202C222.61,178.123,210.721,198.581,216.953,217.982z M266.935,217.404 c-5.501,9.528-17.732,12.802-27.261,7.302c-17.682-10.23-10.301-37.247,10.032-37.247 C264.984,187.459,274.602,204.112,266.935,217.404z"/>
        </g>
    </svg>;
}

export function RulerSvg() {
    return <svg viewBox="0 0 508.073 508.073" width="17px" height="17px">
        <g>
            <path d="M470.459,378.925c-0.7-2.1-1.9-4-3.4-5.5l-113.9-113.9l149.8-149.8c10-10,2.6-17.3,0-19.9l-85.7-85.7 c-5.5-5.5-14.4-5.5-19.9,0l-149.8,149.8l-134.7-134.7c-25.5-25.5-67.8-25.6-93.4,0c-25.8,25.8-25.8,67.7,0,93.5l134.6,134.7 l-149.9,149.9c-5.5,5.5-5.5,14.4,0,19.9l85.6,85.7c2.6,2.6,10,10,19.9,0l150-149.9l113.9,113.9c1.5,1.5,3.4,2.7,5.5,3.4 l110.4,36.9c6,2,10.7,0.3,14.4-3.4c3.8-3.8,5.1-9.4,3.4-14.4L470.459,378.925z M276.159,165.225l29.9,29.9 c5.5,5.5,14.4,5.5,19.9,0c5-5,5.5-14.4,0-19.9l-29.9-29.9l23.7-23.7l13,13c5.5,5.5,14.4,5.5,19.9,0c5.2-5.2,6.2-13.7,0-19.9 l-13-13l23.7-23.7l29.4,29.4c5.5,5.5,14.4,5.5,19.9,0c5.6-5.6,5.5-14.4,0-19.9l-29.4-29.5l24-24l65.8,65.7l-139.8,139.9 l-65.8-65.8L276.159,165.225z M39.359,92.825c-14.8-14.8-14.8-38.8,0-53.6c14.1-14.1,38.8-14.7,53.6,0l15.5,15.5l-53.6,53.6 L39.359,92.825z M99.759,473.025l-65.7-65.7l24-24l13.2,13.2c5.5,5.5,14.4,5.6,19.9,0c5.5-5.5,5.5-14.4,0-19.9l-13.1-13.3 l23.7-23.7l29.6,29.6c6.2,6.2,14.7,5.2,19.9,0c5.5-5.5,5.5-14.4,0-19.9l-29.6-29.6l23.7-23.7l13.2,13.2c5.5,5.5,14.9,5,19.9,0 c5.5-5.5,5.5-14.4,0-19.9l-13.2-13.2l8.8-8.8l65.8,65.7L99.759,473.025z M74.759,128.225l53.6-53.6l308.7,308.7l-53.6,53.6 L74.759,128.225z M409.659,450.725l41.3-41.3l20.7,61.9L409.659,450.725z"/>
        </g>
    </svg>;
}

export function CubeSvg() {
    return <svg viewBox="0 0 270.06 270.06" width="17px" height="17px">
        <g>
            <path d="M264.898,0.007 c-0.181,0.006-0.362,0.023-0.541,0.049H84.996c-1.326,0-2.598,0.527-3.535,1.465l-80,80C0.525,82.459-0.001,83.73,0,85.056v180 c0,2.761,2.239,5,5,5h180c1.326,0,2.598-0.527,3.535-1.465l80-80c0.938-0.938,1.465-2.209,1.465-3.535V5.819 c0.14-0.893,0.035-1.807-0.303-2.645c-0.016-0.045-0.033-0.09-0.051-0.135c-0.808-1.89-2.69-3.093-4.744-3.033L264.898,0.007z M87.066,10.056h165.859l-70,70H17.066L87.066,10.056z M259.996,17.126v165.859l-70,70V87.126L259.996,17.126L259.996,17.126z M84.92,19.985c-2.759,0.042-4.963,2.311-4.924,5.07v40c-0.039,2.761,2.168,5.032,4.929,5.071s5.032-2.168,5.071-4.929 c0.001-0.047,0.001-0.094,0-0.141v-40c0.039-2.761-2.168-5.031-4.93-5.07C85.018,19.984,84.969,19.984,84.92,19.985z M9.996,90.056 h170v170h-170V90.056L9.996,90.056z M84.92,99.985c-2.759,0.042-4.963,2.311-4.924,5.07v30c-0.039,2.761,2.168,5.032,4.929,5.071 s5.032-2.168,5.071-4.929c0.001-0.047,0.001-0.094,0-0.141v-30c0.039-2.761-2.168-5.031-4.93-5.07 C85.018,99.984,84.969,99.984,84.92,99.985z M84.92,159.985c-2.759,0.042-4.963,2.311-4.924,5.07v17.93L61.461,201.52 c-1.992,1.913-2.057,5.078-0.144,7.07c1.913,1.992,5.078,2.057,7.07,0.144c0.049-0.047,0.097-0.095,0.144-0.144l18.535-18.535h17.93 c2.761,0.039,5.032-2.168,5.071-4.929s-2.168-5.032-4.929-5.071c-0.047-0.001-0.094-0.001-0.141,0h-15v-15 c0.039-2.761-2.168-5.031-4.93-5.07C85.018,159.984,84.969,159.984,84.92,159.985z M134.996,180.055 c-2.761-0.039-5.032,2.168-5.071,4.929c-0.039,2.761,2.168,5.032,4.929,5.071c0.047,0.001,0.094,0.001,0.141,0h30 c2.761,0.039,5.032-2.168,5.071-4.929c0.039-2.761-2.168-5.032-4.929-5.071c-0.047-0.001-0.094-0.001-0.141,0H134.996z M204.996,180.055c-2.761-0.039-5.032,2.168-5.071,4.929c-0.039,2.761,2.168,5.032,4.929,5.071c0.047,0.001,0.094,0.001,0.141,0h40 c2.761,0.039,5.032-2.168,5.071-4.929s-2.168-5.032-4.929-5.071c-0.047-0.001-0.094-0.001-0.141,0H204.996z M44.898,220.007 c-1.299,0.039-2.532,0.582-3.438,1.514l-20,20c-1.992,1.913-2.057,5.078-0.144,7.07c1.913,1.992,5.078,2.057,7.07,0.144 c0.049-0.047,0.097-0.095,0.144-0.144l20-20c1.98-1.925,2.025-5.091,0.1-7.071C47.654,220.514,46.3,219.965,44.898,220.007z"/>
        </g>
    </svg>;
}

// The following icons are adapted from https://materialdesignicons.com/ and
// licensed with https://github.com/Templarian/MaterialDesign/blob/master/LICENSE

export function CursorDefaultOutlineSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M10.07,14.27C10.57,14.03 11.16,14.25 11.4,14.75L13.7,19.74L15.5,18.89L13.19,13.91C12.95,13.41 13.17,12.81 13.67,12.58L13.95,12.5L16.25,12.05L8,5.12V15.9L9.82,14.43L10.07,14.27M13.64,21.97C13.14,22.21 12.54,22 12.31,21.5L10.13,16.76L7.62,18.78C7.45,18.92 7.24,19 7,19A1,1 0 0,1 6,18V3A1,1 0 0,1 7,2C7.24,2 7.47,2.09 7.64,2.23L7.65,2.22L19.14,11.86C19.57,12.22 19.62,12.85 19.27,13.27C19.12,13.45 18.91,13.57 18.7,13.61L15.54,14.23L17.74,18.96C18,19.46 17.76,20.05 17.26,20.28L13.64,21.97Z" />
    </svg>;
}

export function FileOutlineSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24" strokeWidth="0.1">
        <path fill="currentColor" d="M14,2H6A2,2 0 0,0 4,4V20A2,2 0 0,0 6,22H18A2,2 0 0,0 20,20V8L14,2M18,20H6V4H13V9H18V20Z" />
    </svg>;
}

// The following icons are adapted from https://material-ui.com/components/material-icons/ and
// licensed with https://github.com/mui-org/material-ui/blob/master/LICENSE

export function AccountTreeOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M22 11V3h-7v3H9V3H2v8h7V8h2v10h4v3h7v-8h-7v3h-2V8h2v3h7zM7 9H4V5h3v4zm10 6h3v4h-3v-4zm0-10h3v4h-3V5z" />
    </svg>;
}

export function AddSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19 13h-6v6h-2v-6H5v-2h6V5h2v6h6v2z" />
    </svg>;
}

export function ArrowDownwardSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M20 12l-1.41-1.41L13 16.17V4h-2v12.17l-5.58-5.59L4 12l8 8 8-8z" />
    </svg>;
}

export function ArrowDropDownSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M7 10l5 5 5-5z" />
    </svg>;
}

export function ArrowRightSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M10 17l5-5-5-5v10z" />
    </svg>;
}

export function ArrowUpwardSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M4 12l1.41 1.41L11 7.83V20h2V7.83l5.58 5.59L20 12l-8-8-8 8z" />
    </svg>;
}

export function AutorenewSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 6v3l4-4-4-4v3c-4.42 0-8 3.58-8 8 0 1.57.46 3.03 1.24 4.26L6.7 14.8c-.45-.83-.7-1.79-.7-2.8 0-3.31 2.69-6 6-6zm6.76 1.74L17.3 9.2c.44.84.7 1.79.7 2.8 0 3.31-2.69 6-6 6v-3l-4 4 4 4v-3c4.42 0 8-3.58 8-8 0-1.57-.46-3.03-1.24-4.26z" />
    </svg>;
}

export function BlurOnSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M6 13c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm0 4c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm0-8c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm-3 .5c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zM6 5c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm15 5.5c.28 0 .5-.22.5-.5s-.22-.5-.5-.5-.5.22-.5.5.22.5.5.5zM14 7c.55 0 1-.45 1-1s-.45-1-1-1-1 .45-1 1 .45 1 1 1zm0-3.5c.28 0 .5-.22.5-.5s-.22-.5-.5-.5-.5.22-.5.5.22.5.5.5zm-11 10c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zm7 7c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zm0-17c.28 0 .5-.22.5-.5s-.22-.5-.5-.5-.5.22-.5.5.22.5.5.5zM10 7c.55 0 1-.45 1-1s-.45-1-1-1-1 .45-1 1 .45 1 1 1zm0 5.5c-.83 0-1.5.67-1.5 1.5s.67 1.5 1.5 1.5 1.5-.67 1.5-1.5-.67-1.5-1.5-1.5zm8 .5c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm0 4c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm0-8c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm0-4c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm3 8.5c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zM14 17c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm0 3.5c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5zm-4-12c-.83 0-1.5.67-1.5 1.5s.67 1.5 1.5 1.5 1.5-.67 1.5-1.5-.67-1.5-1.5-1.5zm0 8.5c-.55 0-1 .45-1 1s.45 1 1 1 1-.45 1-1-.45-1-1-1zm4-4.5c-.83 0-1.5.67-1.5 1.5s.67 1.5 1.5 1.5 1.5-.67 1.5-1.5-.67-1.5-1.5-1.5zm0-4c-.83 0-1.5.67-1.5 1.5s.67 1.5 1.5 1.5 1.5-.67 1.5-1.5-.67-1.5-1.5-1.5z" />
    </svg>;
}

export function BookmarksOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M15 7v12.97l-4.21-1.81-.79-.34-.79.34L5 19.97V7h10m4-6H8.99C7.89 1 7 1.9 7 3h10c1.1 0 2 .9 2 2v13l2 1V3c0-1.1-.9-2-2-2zm-4 4H5c-1.1 0-2 .9-2 2v16l7-3 7 3V7c0-1.1-.9-2-2-2z" />
    </svg>;
}

export function BrushSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M7 14c-1.66 0-3 1.34-3 3 0 1.31-1.16 2-2 2 .92 1.22 2.49 2 4 2 2.21 0 4-1.79 4-4 0-1.66-1.34-3-3-3zm13.71-9.37l-1.34-1.34a.9959.9959 0 00-1.41 0L9 12.25 11.75 15l8.96-8.96c.39-.39.39-1.02 0-1.41z" />
    </svg>;
}

export function BuildOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M22.61 18.99l-9.08-9.08c.93-2.34.45-5.1-1.44-7C9.79.61 6.21.4 3.66 2.26L7.5 6.11 6.08 7.52 2.25 3.69C.39 6.23.6 9.82 2.9 12.11c1.86 1.86 4.57 2.35 6.89 1.48l9.11 9.11c.39.39 1.02.39 1.41 0l2.3-2.3c.4-.38.4-1.01 0-1.41zm-3 1.6l-9.46-9.46c-.61.45-1.29.72-2 .82-1.36.2-2.79-.21-3.83-1.25C3.37 9.76 2.93 8.5 3 7.26l3.09 3.09 4.24-4.24-3.09-3.09c1.24-.07 2.49.37 3.44 1.31 1.08 1.08 1.49 2.57 1.24 3.96-.12.71-.42 1.37-.88 1.96l9.45 9.45-.88.89z" />
    </svg>;
}

export function BuildSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M22.7 19l-9.1-9.1c.9-2.3.4-5-1.5-6.9-2-2-5-2.4-7.4-1.3L9 6 6 9 1.6 4.7C.4 7.1.9 10.1 2.9 12.1c1.9 1.9 4.6 2.4 6.9 1.5l9.1 9.1c.4.4 1 .4 1.4 0l2.3-2.3c.5-.4.5-1.1.1-1.4z" />
    </svg>;
}

export function CameraOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M14.25 2.26l-.08-.04-.01.02C13.46 2.09 12.74 2 12 2 6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10c0-4.75-3.31-8.72-7.75-9.74zM19.41 9h-7.99l2.71-4.7c2.4.66 4.35 2.42 5.28 4.7zM13.1 4.08L10.27 9l-1.15 2L6.4 6.3C7.84 4.88 9.82 4 12 4c.37 0 .74.03 1.1.08zM5.7 7.09L8.54 12l1.15 2H4.26C4.1 13.36 4 12.69 4 12c0-1.85.64-3.55 1.7-4.91zM4.59 15h7.98l-2.71 4.7c-2.4-.67-4.34-2.42-5.27-4.7zm6.31 4.91L14.89 13l2.72 4.7C16.16 19.12 14.18 20 12 20c-.38 0-.74-.04-1.1-.09zm7.4-3l-4-6.91h5.43c.17.64.27 1.31.27 2 0 1.85-.64 3.55-1.7 4.91z" />
    </svg>;
}

export function CameraSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M9.4 10.5l4.77-8.26C13.47 2.09 12.75 2 12 2c-2.4 0-4.6.85-6.32 2.25l3.66 6.35.06-.1zM21.54 9c-.92-2.92-3.15-5.26-6-6.34L11.88 9h9.66zm.26 1h-7.49l.29.5 4.76 8.25C21 16.97 22 14.61 22 12c0-.69-.07-1.35-.2-2zM8.54 12l-3.9-6.75C3.01 7.03 2 9.39 2 12c0 .69.07 1.35.2 2h7.49l-1.15-2zm-6.08 3c.92 2.92 3.15 5.26 6 6.34L12.12 15H2.46zm11.27 0l-3.9 6.76c.7.15 1.42.24 2.17.24 2.4 0 4.6-.85 6.32-2.25l-3.66-6.35-.93 1.6z" />
    </svg>;
}

export function CancelOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 2C6.47 2 2 6.47 2 12s4.47 10 10 10 10-4.47 10-10S17.53 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8zm3.59-13L12 10.59 8.41 7 7 8.41 10.59 12 7 15.59 8.41 17 12 13.41 15.59 17 17 15.59 13.41 12 17 8.41z" />
    </svg>;
}

export function CancelSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 2C6.47 2 2 6.47 2 12s4.47 10 10 10 10-4.47 10-10S17.53 2 12 2zm5 13.59L15.59 17 12 13.41 8.41 17 7 15.59 10.59 12 7 8.41 8.41 7 12 10.59 15.59 7 17 8.41 13.41 12 17 15.59z" />
    </svg>;
}

export function CenterFocusStrongSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 8c-2.21 0-4 1.79-4 4s1.79 4 4 4 4-1.79 4-4-1.79-4-4-4zm-7 7H3v4c0 1.1.9 2 2 2h4v-2H5v-4zM5 5h4V3H5c-1.1 0-2 .9-2 2v4h2V5zm14-2h-4v2h4v4h2V5c0-1.1-.9-2-2-2zm0 16h-4v2h4c1.1 0 2-.9 2-2v-4h-2v4z" />
    </svg>;
}

export function CheckSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z" />
    </svg>;
}

export function ClearSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z" />
    </svg>;
}

export function CloseSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z" />
    </svg>;
}

export function CloudUploadSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19.35 10.04C18.67 6.59 15.64 4 12 4 9.11 4 6.6 5.64 5.35 8.04 2.34 8.36 0 10.91 0 14c0 3.31 2.69 6 6 6h13c2.76 0 5-2.24 5-5 0-2.64-2.05-4.78-4.65-4.96zM14 13v4h-4v-4H7l5-5 5 5h-3z" />
    </svg>;
}

export function CodeSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M9.4 16.6L4.8 12l4.6-4.6L8 6l-6 6 6 6 1.4-1.4zm5.2 0l4.6-4.6-4.6-4.6L16 6l6 6-6 6-1.4-1.4z" />
    </svg>;
}

export function DeleteOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M16 9v10H8V9h8m-1.5-6h-5l-1 1H5v2h14V4h-3.5l-1-1zM18 7H6v12c0 1.1.9 2 2 2h8c1.1 0 2-.9 2-2V7z" />
    </svg>;
}

export function DeleteSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M6 19c0 1.1.9 2 2 2h8c1.1 0 2-.9 2-2V7H6v12zM19 4h-3.5l-1-1h-5l-1 1H5v2h14V4z" />
    </svg>;
}

export function ErrorSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z" />
    </svg>;
}

export function ExtensionSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M20.5 11H19V7c0-1.1-.9-2-2-2h-4V3.5C13 2.12 11.88 1 10.5 1S8 2.12 8 3.5V5H4c-1.1 0-1.99.9-1.99 2v3.8H3.5c1.49 0 2.7 1.21 2.7 2.7s-1.21 2.7-2.7 2.7H2V20c0 1.1.9 2 2 2h3.8v-1.5c0-1.49 1.21-2.7 2.7-2.7 1.49 0 2.7 1.21 2.7 2.7V22H17c1.1 0 2-.9 2-2v-4h1.5c1.38 0 2.5-1.12 2.5-2.5S21.88 11 20.5 11z" />
    </svg>;
}

export function FlipToFrontSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24" strokeWidth="0.1px">
        <path d="M3 13h2v-2H3v2zm0 4h2v-2H3v2zm2 4v-2H3c0 1.1.89 2 2 2zM3 9h2V7H3v2zm12 12h2v-2h-2v2zm4-18H9c-1.11 0-2 .9-2 2v10c0 1.1.89 2 2 2h10c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm0 12H9V5h10v10zm-8 6h2v-2h-2v2zm-4 0h2v-2H7v2z" />
    </svg>;
}

export function FullscreenSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M7 14H5v5h5v-2H7v-3zm-2-4h2V7h3V5H5v5zm12 7h-3v2h5v-5h-2v3zM14 5v2h3v3h2V5h-5z" />
    </svg>;
}

export function GetAppSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z" />
    </svg>;
}

export function HelpOutlineSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M11 18h2v-2h-2v2zm1-16C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8zm0-14c-2.21 0-4 1.79-4 4h2c0-1.1.9-2 2-2s2 .9 2 2c0 2-3 1.75-3 5h2c0-2.25 3-2.5 3-5 0-2.21-1.79-4-4-4z" />
    </svg>;
}

export function HomeOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 5.69l5 4.5V18h-2v-6H9v6H7v-7.81l5-4.5M12 3L2 12h3v8h6v-6h2v6h6v-8h3L12 3z" />
    </svg>;
}

export function LaunchSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19 19H5V5h7V3H5c-1.11 0-2 .9-2 2v14c0 1.1.89 2 2 2h14c1.1 0 2-.9 2-2v-7h-2v7zM14 3v2h3.59l-9.83 9.83 1.41 1.41L19 6.41V10h2V3h-7z" />
    </svg>;
}

export function LinearScaleSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19.5 9.5c-1.03 0-1.9.62-2.29 1.5h-2.92c-.39-.88-1.26-1.5-2.29-1.5s-1.9.62-2.29 1.5H6.79c-.39-.88-1.26-1.5-2.29-1.5C3.12 9.5 2 10.62 2 12s1.12 2.5 2.5 2.5c1.03 0 1.9-.62 2.29-1.5h2.92c.39.88 1.26 1.5 2.29 1.5s1.9-.62 2.29-1.5h2.92c.39.88 1.26 1.5 2.29 1.5 1.38 0 2.5-1.12 2.5-2.5s-1.12-2.5-2.5-2.5z" />
    </svg>;
}

export function MoreHorizSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M6 10c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm12 0c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm-6 0c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2z" />
    </svg>;
}

export function NavigateBeforeSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M15.41 7.41L14 6l-6 6 6 6 1.41-1.41L10.83 12z" />
    </svg>;
}

export function NavigateNextSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M10 6L8.59 7.41 13.17 12l-4.58 4.59L10 18l6-6z" />
    </svg>;
}

export function OpenInBrowserSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19 4H5c-1.11 0-2 .9-2 2v12c0 1.1.89 2 2 2h4v-2H5V8h14v10h-4v2h4c1.1 0 2-.9 2-2V6c0-1.1-.89-2-2-2zm-7 6l-4 4h3v6h2v-6h3l-4-4z" />
    </svg>;
}

export function PlayArrowSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M8 5v14l11-7z" />
    </svg>;
}

export function RefreshSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M17.65 6.35C16.2 4.9 14.21 4 12 4c-4.42 0-7.99 3.58-7.99 8s3.57 8 7.99 8c3.73 0 6.84-2.55 7.73-6h-2.08c-.82 2.33-3.04 4-5.65 4-3.31 0-6-2.69-6-6s2.69-6 6-6c1.66 0 3.14.69 4.22 1.78L13 11h7V4l-2.35 2.35z" />
    </svg>;
}

export function RemoveSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M19 13H5v-2h14v2z" />
    </svg>;
}

export function RestoreSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M13 3c-4.97 0-9 4.03-9 9H1l3.89 3.89.07.14L9 12H6c0-3.87 3.13-7 7-7s7 3.13 7 7-3.13 7-7 7c-1.93 0-3.68-.79-4.94-2.06l-1.42 1.42C8.27 19.99 10.51 21 13 21c4.97 0 9-4.03 9-9s-4.03-9-9-9zm-1 5v5l4.28 2.54.72-1.21-3.5-2.08V8H12z" />
    </svg>;
}

export function SaveOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M17 3H5c-1.11 0-2 .9-2 2v14c0 1.1.89 2 2 2h14c1.1 0 2-.9 2-2V7l-4-4zm2 16H5V5h11.17L19 7.83V19zm-7-7c-1.66 0-3 1.34-3 3s1.34 3 3 3 3-1.34 3-3-1.34-3-3-3zM6 6h9v4H6z" />
    </svg>;
}

export function ScatterPlotSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <circle cx="7" cy="14" r="3" /><circle cx="11" cy="6" r="3" /><circle cx="16.6" cy="17.6" r="3" />
    </svg>;
}

export function SkipPreviousSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M6 6h2v12H6zm3.5 6l8.5 6V6z" />
    </svg>;
}

export function StopSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M6 6h12v12H6z" />
    </svg>;
}

export function SubscriptionsOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M4 6h16v2H4zm2-4h12v2H6zm14 8H4c-1.1 0-2 .9-2 2v8c0 1.1.9 2 2 2h16c1.1 0 2-.9 2-2v-8c0-1.1-.9-2-2-2zm0 10H4v-8h16v8zm-10-7.27v6.53L16 16z" />
    </svg>;
}

export function SwapHorizSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M6.99 11L3 15l3.99 4v-3H14v-2H6.99v-3zM21 9l-3.99-4v3H10v2h7.01v3L21 9z" />
    </svg>;
}

export function TuneSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M3 17v2h6v-2H3zM3 5v2h10V5H3zm10 16v-2h8v-2h-8v-2h-2v6h2zM7 9v2H3v2h4v2h2V9H7zm14 4v-2H11v2h10zm-6-4h2V7h4V5h-4V3h-2v6z" />
    </svg>;
}

export function VisibilityOffOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 6c3.79 0 7.17 2.13 8.82 5.5-.59 1.22-1.42 2.27-2.41 3.12l1.41 1.41c1.39-1.23 2.49-2.77 3.18-4.53C21.27 7.11 17 4 12 4c-1.27 0-2.49.2-3.64.57l1.65 1.65C10.66 6.09 11.32 6 12 6zm-1.07 1.14L13 9.21c.57.25 1.03.71 1.28 1.28l2.07 2.07c.08-.34.14-.7.14-1.07C16.5 9.01 14.48 7 12 7c-.37 0-.72.05-1.07.14zM2.01 3.87l2.68 2.68C3.06 7.83 1.77 9.53 1 11.5 2.73 15.89 7 19 12 19c1.52 0 2.98-.29 4.32-.82l3.42 3.42 1.41-1.41L3.42 2.45 2.01 3.87zm7.5 7.5l2.61 2.61c-.04.01-.08.02-.12.02-1.38 0-2.5-1.12-2.5-2.5 0-.05.01-.08.01-.13zm-3.4-3.4l1.75 1.75c-.23.55-.36 1.15-.36 1.78 0 2.48 2.02 4.5 4.5 4.5.63 0 1.23-.13 1.77-.36l.98.98c-.88.24-1.8.38-2.75.38-3.79 0-7.17-2.13-8.82-5.5.7-1.43 1.72-2.61 2.93-3.53z" />
    </svg>;
}

export function VisibilityOutlinedSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M12 6c3.79 0 7.17 2.13 8.82 5.5C19.17 14.87 15.79 17 12 17s-7.17-2.13-8.82-5.5C4.83 8.13 8.21 6 12 6m0-2C7 4 2.73 7.11 1 11.5 2.73 15.89 7 19 12 19s9.27-3.11 11-7.5C21.27 7.11 17 4 12 4zm0 5c1.38 0 2.5 1.12 2.5 2.5S13.38 14 12 14s-2.5-1.12-2.5-2.5S10.62 9 12 9m0-2c-2.48 0-4.5 2.02-4.5 4.5S9.52 16 12 16s4.5-2.02 4.5-4.5S14.48 7 12 7z" />
    </svg>;
}

export function WarningSvg() {
    return <svg width="24px" height="24px" viewBox="0 0 24 24">
        <path d="M1 21h22L12 2 1 21zm12-3h-2v-2h2v2zm0-4h-2v-4h2v4z" />
    </svg>;
}

// Aliases

export const SelectionModeSvg = CursorDefaultOutlineSvg;
export const SuperposeAtomsSvg = ScatterPlotSvg;
export const SuperposeChainsSvg = LinearScaleSvg;
export const SuperpositionSvg = FlipToFrontSvg;
