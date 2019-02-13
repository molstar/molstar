/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from './context';
import { Plugin } from './ui/plugin'
import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { PluginCommands } from './command';
import { PluginSpec } from './spec';
import { DownloadStructure, CreateComplexRepresentation, OpenStructure, OpenVolume, DownloadDensity } from './state/actions/basic';
import { StateTransforms } from './state/transforms';
import { PluginBehaviors } from './behavior';

function getParam(name: string, regex: string): string {
    let r = new RegExp(`${name}=(${regex})[&]?`, 'i');
    return decodeURIComponent(((window.location.search || '').match(r) || [])[1] || '');
}

export const DefaultPluginSpec: PluginSpec = {
    actions: [
        PluginSpec.Action(DownloadStructure),
        PluginSpec.Action(DownloadDensity),
        PluginSpec.Action(OpenStructure),
        PluginSpec.Action(OpenVolume),
        PluginSpec.Action(CreateComplexRepresentation),
        PluginSpec.Action(StateTransforms.Data.Download),
        PluginSpec.Action(StateTransforms.Data.ParseCif),
        PluginSpec.Action(StateTransforms.Data.ParseCcp4),
        PluginSpec.Action(StateTransforms.Model.StructureAssemblyFromModel),
        PluginSpec.Action(StateTransforms.Model.StructureSymmetryFromModel),
        PluginSpec.Action(StateTransforms.Model.StructureFromModel),
        PluginSpec.Action(StateTransforms.Model.ModelFromTrajectory),
        PluginSpec.Action(StateTransforms.Model.VolumeFromCcp4),
        PluginSpec.Action(StateTransforms.Representation.StructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.ExplodeStructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.VolumeRepresentation3D),
    ],
    behaviors: [
        PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
        PluginSpec.Behavior(PluginBehaviors.Camera.FocusLociOnSelect, { minRadius: 20, extraRadius: 4 }),
        PluginSpec.Behavior(PluginBehaviors.Animation.StructureAnimation, { rotate: false, rotateValue: 0, explode: false, explodeValue: 0 }),
        PluginSpec.Behavior(PluginBehaviors.Labels.SceneLabels),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.PDBeStructureQualityReport, { autoAttach: true }),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.RCSBAssemblySymmetry, { autoAttach: true }),
    ]
}

export function createPlugin(target: HTMLElement, spec?: PluginSpec): PluginContext {
    const ctx = new PluginContext(spec || DefaultPluginSpec);
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
        ctx.log.error('Failed to load snapshot.');
        console.warn('Failed to load snapshot', e);
    }
}