/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { Plugin } from './plugin';
import { PluginUIContext } from './context';
import { DefaultPluginUISpec, PluginUISpec } from './spec';

export async function createPluginUI(target: HTMLElement, spec?: PluginUISpec, options?: { onBeforeUIRender?: (ctx: PluginUIContext) => (Promise<void> | void) }) {
    const ctx = new PluginUIContext(spec || DefaultPluginUISpec());
    await ctx.init();
    if (options?.onBeforeUIRender) {
        await options.onBeforeUIRender(ctx);
    }
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    return ctx;
}