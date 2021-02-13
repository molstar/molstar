/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { Plugin } from '../mol-plugin-ui/plugin';
import { PluginContext } from './context';
import { DefaultPluginSpec, PluginSpec } from './spec';

export function createPlugin(target: HTMLElement, spec?: PluginSpec): PluginContext {
    const ctx = new PluginContext(spec || DefaultPluginSpec());
    ctx.init();
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    return ctx;
}

/** Returns the instance of the plugin after all behaviors have been initialized */
export async function createPluginAsync(target: HTMLElement, spec?: PluginSpec) {
    const ctx = new PluginContext(spec || DefaultPluginSpec());
    const init = ctx.init();
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    await init;
    return ctx;
}