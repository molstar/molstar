/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createPlugin, DefaultPluginSpec } from '../../mol-plugin';
import './index.html'
import { PluginContext } from '../../mol-plugin/context';
import { PluginCommands } from '../../mol-plugin/command';
import { PluginSpec } from '../../mol-plugin/spec';
import { CreateJoleculeState } from './extensions/jolecule';
require('mol-plugin/skin/light.scss')

function getParam(name: string, regex: string): string {
    let r = new RegExp(`${name}=(${regex})[&]?`, 'i');
    return decodeURIComponent(((window.location.search || '').match(r) || [])[1] || '');
}

const hideControls = getParam('hide-controls', `[^&]+`) === '1';

function init() {
    const spec: PluginSpec = {
        actions: [...DefaultPluginSpec.actions, PluginSpec.Action(CreateJoleculeState)],
        behaviors: [...DefaultPluginSpec.behaviors],
        animations: [...DefaultPluginSpec.animations || []],
        customParamEditors: DefaultPluginSpec.customParamEditors,
        layout: {
            initial: {
                isExpanded: true,
                showControls: !hideControls
            },
            controls: {
                ...DefaultPluginSpec.layout && DefaultPluginSpec.layout.controls
            }
        }
    };
    const plugin = createPlugin(document.getElementById('app')!, spec);
    trySetSnapshot(plugin);
}

async function trySetSnapshot(ctx: PluginContext) {
    try {
        const snapshotUrl = getParam('snapshot-url', `[^&]+`);
        const snapshotId = getParam('snapshot-id', `[^&]+`);
        if (!snapshotUrl && !snapshotId) return;
        // TODO parametrize the server
        const url = snapshotId
            ? `https://webchem.ncbr.muni.cz/molstar-state/get/${snapshotId}`
            : snapshotUrl;
        await PluginCommands.State.Snapshots.Fetch.dispatch(ctx, { url })
    } catch (e) {
        ctx.log.error('Failed to load snapshot.');
        console.warn('Failed to load snapshot', e);
    }
}

init();