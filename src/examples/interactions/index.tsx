/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComputeInteractions } from '../../extensions/interactions/transforms';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import '../../mol-plugin-ui/skin/light.scss';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginSpec } from '../../mol-plugin/spec';
import './index.html';

async function createViewer(root: HTMLElement) {
    const spec = DefaultPluginUISpec();
    const plugin = await createPluginUI({
        target: root,
        render: renderReact18,
        spec: {
            ...spec,
            layout: {
                initial: {
                    isExpanded: true,
                    showControls: false
                }
            },
            components: {
                remoteState: 'none',
            },
            behaviors: [
                ...spec.behaviors,
                PluginSpec.Behavior(MolViewSpec)
            ],
            config: [
                [PluginConfig.Viewport.ShowAnimation, false],
                [PluginConfig.Viewport.ShowTrajectoryControls, false],
            ]
        }
    });

    return plugin;
}

export async function loadInteractionsExample(root: HTMLElement) {
    const plugin = await createViewer(root);

    const receptorData = await plugin.builders.data.download({ url: '../../../examples/ace2.pdbqt' });
    const receptorTrajectory = await plugin.builders.structure.parseTrajectory(receptorData, 'pdbqt');
    const receptor = await plugin.builders.structure.hierarchy.applyPreset(receptorTrajectory, 'default');

    const ligandData = await plugin.builders.data.download({ url: '../../../examples/ace2-hit.mol2' });
    const ligandTrajectory = await plugin.builders.structure.parseTrajectory(ligandData, 'mol2');
    const ligand = await plugin.builders.structure.hierarchy.applyPreset(ligandTrajectory, 'default', { representationPreset: 'atomic-detail' });

    const update = plugin.build();

    const refs = [
        receptor?.representation.components.polymer.ref,
        ligand?.representation.components.all.ref,
    ];

    const interactions = update.toRoot().apply(ComputeInteractions, {
        sources: refs.map(structureRef => ({ structureRef })),
    }, { dependsOn: refs }).selector;

    await update.commit();

    console.log(interactions.data);
}

loadInteractionsExample(document.getElementById('viewer')!);
