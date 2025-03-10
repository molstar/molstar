/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureInteractions } from '../../extensions/interactions/model';
import { ComputeInteractions, InteractionsShape } from '../../extensions/interactions/transforms';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { forEachSchemaItem, StructureElementSchema, StructureElementSchemaItem, structureElementSchemaToExpression } from '../../mol-model/structure/query/schema';
import { StructureSelectionQuery } from '../../mol-plugin-state/helpers/structure-selection-query';
import { ShapeRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import '../../mol-plugin-ui/skin/light.scss';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginCommands } from '../../mol-plugin/commands';
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

    const receptorRef = receptor?.representation.components.polymer.ref!;
    const ligandRef = ligand?.representation.components.all.ref!;

    const refs = [receptorRef, ligandRef];

    const interactionsRef = update.toRoot().apply(ComputeInteractions, {
        sources: refs.map(structureRef => ({ structureRef })),
    }, { dependsOn: refs });

    interactionsRef.apply(InteractionsShape).apply(ShapeRepresentation3D);

    await update.commit();

    const expression = getInteractionExpression(interactionsRef.selector.data?.interactions!, receptorRef);
    await plugin.managers.structure.component.add({
        selection: StructureSelectionQuery('Binding Site', expression),
        representation: 'ball-and-stick',
        options: { checkExisting: false, label: 'Binding Site' },
    });

    PluginCommands.Camera.FocusObject(plugin, {
        targets: [{
            targetRef: ligand?.representation.representations.all.ref
        }]
    });
}

function getInteractionExpression(interactions: StructureInteractions, receptorRef: string) {
    // TODO: check if empty, etc.

    const added = new Set<string>();
    const schema: StructureElementSchemaItem[] = [];

    const add = (a: StructureElementSchema) => {
        forEachSchemaItem(a, item => {
            const key = JSON.stringify(item);
            if (!added.has(key)) {
                added.add(key);
                schema.push(item);
            }
        });
    };

    for (const e of interactions.elements) {
        // TODO: add the residues from the LOCI instead!
        if (e.schema.aStructureRef === receptorRef) add(normalizeSchema(e.schema.a));
        if (e.schema.bStructureRef === receptorRef) add(normalizeSchema(e.schema.b));
    }

    return structureElementSchemaToExpression(schema);
}

function normalizeSchema(schema: StructureElementSchema) {
    const ret: StructureElementSchema = [];
    forEachSchemaItem(schema, item => {
        const e = { ...item };
        delete e.atom_id;
        delete e.atom_index;
        delete e.label_atom_id;
        delete e.auth_atom_id;
        ret.push(e);
    });
    return ret;
}

loadInteractionsExample(document.getElementById('viewer')!);
