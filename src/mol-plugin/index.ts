/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from './context';
import { Plugin } from './ui/plugin'
import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { PluginCommands } from './command';
import { PluginSpec } from './spec';
import { DownloadAtomicStructure, CreateComplexRepresentation, OpenAtomicStructure } from './state/actions/basic';
import { StateTransforms } from './state/transforms';
import { PluginBehaviors } from './behavior';
import { LogEntry } from 'mol-util/log-entry';

function getParam(name: string, regex: string): string {
    let r = new RegExp(`${name}=(${regex})[&]?`, 'i');
    return decodeURIComponent(((window.location.search || '').match(r) || [])[1] || '');
}

const DefaultSpec: PluginSpec = {
    actions: [
        PluginSpec.Action(DownloadAtomicStructure),
        PluginSpec.Action(OpenAtomicStructure),
        PluginSpec.Action(CreateComplexRepresentation),
        PluginSpec.Action(StateTransforms.Data.Download),
        PluginSpec.Action(StateTransforms.Data.ParseCif),
        PluginSpec.Action(StateTransforms.Model.StructureAssemblyFromModel),
        PluginSpec.Action(StateTransforms.Model.StructureFromModel),
        PluginSpec.Action(StateTransforms.Model.ModelFromTrajectory),
        PluginSpec.Action(StateTransforms.Representation.StructureRepresentation3D)
    ],
    behaviors: [
        PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
        PluginSpec.Behavior(PluginBehaviors.Camera.FocusLociOnSelect, { minRadius: 20, extraRadius: 4 })
    ]
}

export function createPlugin(target: HTMLElement): PluginContext {
    const ctx = new PluginContext(DefaultSpec);
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);

    trySetSnapshot(ctx);

    return ctx;
}

async function trySetSnapshot(ctx: PluginContext) {
    try {
        const snapshotUrl = getParam('snapshot-url', `[^&]+`);
        if (!snapshotUrl) return;
        await PluginCommands.State.Snapshots.Fetch.dispatch(ctx, { url: snapshotUrl })
    } catch (e) {
        ctx.log(LogEntry.error('Failed to load snapshot.'));
        console.warn('Failed to load snapshot', e);
    }
}