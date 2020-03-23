/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';

export type IconName = FontIconName | SvgIconName

export function Icon(props: {
    name: IconName | undefined,
    style?: React.CSSProperties,
    title?: string
}) {
    if (!props.name) return null;
    switch (props.name) {
        case 'union':
        case 'subtract':
        case 'intersect':
        case 'set':
            return <SvgIcon name={props.name} title={props.title} style={props.style} />
        default: return <FontIcon name={props.name} title={props.title} style={props.style} />
    }
}

//

type FontIconName =
    | '' | 'expand-layout' | 'plus' | 'minus' | 'reset-scene' | 'ok' | 'back' | 'block' | 'off' | 'expand' | 'collapse' | 'visual-visibility'
    | 'abort' | 'focus-on-visual' | 'settings' | 'tools' | 'log' | 'remove' | 'help' | 'help-circle' | 'info' | 'left-open-big' | 'right-open-big'
    | 'left-open' | 'right-open' | 'screenshot' | 'model-prev' | 'model-next' | 'model-first' | 'down-thin' | 'up-thin' | 'left-thin' | 'right-thin'
    | 'switch' | 'play' | 'stop' | 'pause' | 'cw' | 'database' | 'upload' | 'record' | 'code' | 'floppy' | 'tape' | 'flow-cascade' | 'flow-tree'
    | 'home' | 'address' | 'download' | 'export' | 'palette' | 'search' | 'flashlight' | 'mail' | 'heart' | 'heart-empty' | 'star' | 'star-empty'
    | 'user' | 'users' | 'user-add' | 'video' | 'picture' | 'camera' | 'layout' | 'menu' | 'check' | 'cancel' | 'cancel-circled' | 'cancel-squared'
    | 'plus-circled' | 'plus-squared' | 'minus-circled' | 'minus-squared' | 'help-circled' | 'info-circled' | 'link' | 'attach' | 'lock' | 'lock-open'
    | 'eye' | 'tag' | 'bookmark' | 'bookmarks' | 'flag' | 'thumbs-up' | 'thumbs-down' | 'upload-cloud' | 'reply' | 'reply-all' | 'forward' | 'quote'
    | 'pencil' | 'feather' | 'print' | 'retweet' | 'keyboard' | 'comment' | 'chat' | 'bell' | 'attention' | 'alert' | 'vcard' | 'location' | 'map'
    | 'direction' | 'compass' | 'cup' | 'trash' | 'doc' | 'docs' | 'doc-landscape' | 'doc-text' | 'doc-text-inv' | 'newspaper' | 'book-open' | 'book'
    | 'folder' | 'archive' | 'box' | 'rss' | 'phone' | 'cog' | 'share' | 'shareable' | 'basket' | 'bag' | 'calendar' | 'login' | 'logout' | 'mic'
    | 'mute' | 'sound' | 'volume' | 'clock' | 'hourglass' | 'lamp' | 'light-down' | 'light-up' | 'adjust' | 'resize-full' | 'resize-small' | 'popup'
    | 'publish' | 'window' | 'arrow-combo' | 'down-circled' | 'left-circled' | 'right-circled' | 'up-circled' | 'down-open' | 'up-open' | 'down-open-mini'
    | 'left-open-mini' | 'right-open-mini' | 'up-open-mini' | 'down-open-big' | 'up-open-big' | 'down' | 'left' | 'right' | 'up' | 'down-dir' | 'left-dir'
    | 'right-dir' | 'up-dir' | 'down-bold' | 'left-bold' | 'right-bold' | 'up-bold' | 'ccw' | 'arrows-ccw' | 'level-down' | 'level-up' | 'shuffle'
    | 'loop' | 'to-end' | 'to-start' | 'fast-forward' | 'fast-backward' | 'progress-0' | 'progress-1' | 'progress-2' | 'progress-3' | 'target' | 'list'
    | 'list-add' | 'battery' | 'back-in-time' | 'monitor' | 'mobile' | 'cd' | 'inbox' | 'install' | 'globe' | 'cloud' | 'cloud-thunder' | 'flash'
    | 'moon' | 'flight' | 'paper-plane' | 'leaf' | 'lifebuoy' | 'mouse' | 'briefcase' | 'suitcase' | 'dot' | 'dot-2' | 'dot-3' | 'brush' | 'infinity'
    | 'erase' | 'chart-pie' | 'chart-line' | 'chart-bar' | 'chart-area' | 'graduation-cap' | 'language' | 'ticket' | 'water' | 'droplet' | 'air'
    | 'credit-card' | 'clipboard' | 'megaphone' | 'drive' | 'bucket' | 'thermometer' | 'key' | 'flow-branch' | 'flow-line' | 'flow-parallel' | 'rocket'
    | 'gauge' | 'help-circle-collapse' | 'help-circle-expand'

function FontIcon(props: {
    name: FontIconName,
    style?: React.CSSProperties,
    title?: string
}) {
    return <span className={`msp-icon msp-icon-${props.name}`} style={props.style} title={props.title} />;
}

//

type SvgIconName =
    | '' | 'set' | 'intersect' | 'union' | 'subtract'

function SvgIcon(props: {
    name: SvgIconName,
    style?: React.CSSProperties,
    title?: string
}) {
    return <div className='msp-icon msp-svg-icon' style={props.style} title={props.title}>{getSvg(props.name)}</div>;
}

function getSvg(name: SvgIconName) {
    switch (name) {
        case 'union': return <Union />
        case 'subtract': return <Subtract />
        case 'intersect': return <Intersect />
        case 'set': return <Set />
    }
}

const circleLeft = <circle r="6px" id="circle-left" cy="16px" cx="12px" strokeWidth="1"/>
const circleRight = <circle r="6px" id="circle-right" cy="16px" cx="20px" strokeWidth="1"/>

function Union() {
    return <svg width="32px" height="32px" viewBox="0 0 32 32">
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
    return <svg width="32px" height="32px" viewBox="0 0 32 32">
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
    return <svg width="32px" height="32px" viewBox="0 0 32 32">
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
    return <svg width="32px" height="32px" viewBox="0 0 32 32">
        <defs>
            {circleLeft}
            {circleRight}
        </defs>
        <g>
            <use href="#circle-left" className="msp-shape-empty"/>
            <use href="#circle-right" className="msp-shape-filled"/>
        </g>
    </svg>;
}