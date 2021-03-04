/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { Plugin } from './plugin';
import { PluginUIContext } from './context';
import { DefaultPluginUISpec, PluginUISpec } from './spec';

export function createPlugin(target: HTMLElement, spec?: PluginUISpec): PluginUIContext {
    const ctx = new PluginUIContext(spec || DefaultPluginUISpec());
    ctx.init();
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    return ctx;
}

/** Returns the instance of the plugin after all behaviors have been initialized */
export async function createPluginAsync(target: HTMLElement, spec?: PluginUISpec) {
    const ctx = new PluginUIContext(spec || DefaultPluginUISpec());
    const init = ctx.init();
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    await init;
    return ctx;
}