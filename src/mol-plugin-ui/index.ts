/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */


import { createElement } from 'react';
import { Plugin } from './plugin';
import { PluginUIContext } from './context';
import { DefaultPluginUISpec, PluginUISpec } from './spec';

export async function createPluginUI(options: { target: HTMLElement, render: (component: any, container: Element) => any, spec?: PluginUISpec, onBeforeUIRender?: (ctx: PluginUIContext) => (Promise<void> | void) }) {
    const { spec, target, onBeforeUIRender, render } = options;
    const ctx = new PluginUIContext(spec || DefaultPluginUISpec());
    await ctx.init();
    if (onBeforeUIRender) {
        await onBeforeUIRender(ctx);
    }
    render(createElement(Plugin, { plugin: ctx }), target);
    try {
        await ctx.canvas3dInitialized;
    } catch {
        // Error reported in UI/console elsewhere.
    }
    return ctx;
}