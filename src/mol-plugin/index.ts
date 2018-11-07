/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from './context';
import { Plugin } from './ui/plugin'
import * as React from 'react';
import * as ReactDOM from 'react-dom';

function getParam(name: string, regex: string): string {
    let r = new RegExp(`${name}=(${regex})[&]?`, 'i');
    return decodeURIComponent(((window.location.search || '').match(r) || [])[1] || '');
}

export function createPlugin(target: HTMLElement): PluginContext {
    const ctx = new PluginContext();
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);

    try {
        trySetSnapshot(ctx);
    } catch (e) {
        console.warn('Failed to load snapshot', e);
    }

    return ctx;
}

function trySetSnapshot(ctx: PluginContext) {
    const snapshot = getParam('snapshot', `(?:[A-Za-z0-9+/]{4})*(?:[A-Za-z0-9+/]{2}==|[A-Za-z0-9+/]{3}=)?`);
    if (!snapshot) return;
    const data = JSON.parse(atob(snapshot));
    setTimeout(() => ctx.state.setSnapshot(data), 250);
}