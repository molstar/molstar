/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import '../../mol-util/polyfill';
import { createPlugin, DefaultPluginSpec } from '../../mol-plugin';
import './index.html'
import './favicon.ico'
import { PluginContext } from '../../mol-plugin/context';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginSpec } from '../../mol-plugin/spec';
import { LoadCellPackModel } from './extensions/cellpack/model';
import { StructureFromCellpack } from './extensions/cellpack/state';
import { DownloadStructure } from '../../mol-plugin-state/actions/structure';
require('mol-plugin-ui/skin/light.scss')

function getParam(name: string, regex: string): string {
    let r = new RegExp(`${name}=(${regex})[&]?`, 'i');
    return decodeURIComponent(((window.location.search || '').match(r) || [])[1] || '');
}

const hideControls = getParam('hide-controls', `[^&]+`) === '1';

function init() {
    const spec: PluginSpec = {
        actions: [
            ...DefaultPluginSpec.actions,
            // PluginSpec.Action(CreateJoleculeState),
            PluginSpec.Action(LoadCellPackModel),
            PluginSpec.Action(StructureFromCellpack),
        ],
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
        },
        config: DefaultPluginSpec.config
    };
    const plugin = createPlugin(document.getElementById('app')!, spec);
    trySetSnapshot(plugin);
    tryLoadFromUrl(plugin);
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
        await PluginCommands.State.Snapshots.Fetch(ctx, { url })
    } catch (e) {
        ctx.log.error('Failed to load snapshot.');
        console.warn('Failed to load snapshot', e);
    }
}

async function tryLoadFromUrl(ctx: PluginContext) {
    const url = getParam('loadFromURL', '[^&]+').trim();
    try {
        if (!url) return;

        let format = 'cif', isBinary = false;
        switch (getParam('loadFromURLFormat', '[a-z]+').toLocaleLowerCase().trim()) {
            case 'pdb': format = 'pdb'; break;
            case 'mmbcif': isBinary = true; break;
        }

        const params = DownloadStructure.createDefaultParams(void 0 as any, ctx);

        return ctx.runTask(ctx.state.dataState.applyAction(DownloadStructure, {
            source: {
                name: 'url',
                params: {
                    url,
                    format: format as any,
                    isBinary,
                    options: params.source.params.options,
                    structure: params.source.params.structure,
                }
            }
        }));
    } catch (e) {
        ctx.log.error(`Failed to load from URL (${url})`);
        console.warn(`Failed to load from URL (${url})`, e);
    }
}

init();