/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createPlugin, DefaultPluginSpec } from 'mol-plugin';
import './index.html'
import { PluginContext } from 'mol-plugin/context';
import { PluginCommands } from 'mol-plugin/command';
require('mol-plugin/skin/light.scss')

function getParam(name: string, regex: string): string {
    let r = new RegExp(`${name}=(${regex})[&]?`, 'i');
    return decodeURIComponent(((window.location.search || '').match(r) || [])[1] || '');
}

const hideControls = getParam('hide-controls', `[^&]+`) === '1';

function init() {
    const plugin = createPlugin(document.getElementById('app')!, {
        ...DefaultPluginSpec,
        layout: {
            initial: {
                isExpanded: true,
                showControls: !hideControls
            }
        }
    });
    trySetSnapshot(plugin);
}

async function trySetSnapshot(ctx: PluginContext) {
    try {
        const snapshotUrl = getParam('snapshot-url', `[^&]+`);
        if (!snapshotUrl) return;
        await PluginCommands.State.Snapshots.Fetch.dispatch(ctx, { url: snapshotUrl })
    } catch (e) {
        ctx.log.error('Failed to load snapshot.');
        console.warn('Failed to load snapshot', e);
    }
}

init();