/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from './context';
import { Plugin } from './ui/plugin'
import * as React from 'react';
import * as ReactDOM from 'react-dom';

export function createPlugin(target: HTMLElement): PluginContext {
    const ctx = new PluginContext();
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    return ctx;
}