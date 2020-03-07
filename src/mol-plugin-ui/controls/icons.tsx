/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';

export type IconName = 
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

export function Icon(props: {
    name: IconName | undefined,
    style?: React.CSSProperties,
    title?: string
}) {
    if (!props.name) return null;
    return <span className={`msp-icon msp-icon-${props.name}`} style={props.style} title={props.title} />;
}