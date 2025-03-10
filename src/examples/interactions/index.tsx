/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureInteractions } from '../../extensions/interactions/model';
import { ComputeInteractions, CustomInteractions, InteractionsShape } from '../../extensions/interactions/transforms';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { ResidueIndex, Structure, StructureElement, StructureProperties, StructureQuery } from '../../mol-model/structure';
import { atoms } from '../../mol-model/structure/query/queries/generators';
import { StructureSelectionFromBundle } from '../../mol-plugin-state/transforms/model';
import { ShapeRepresentation3D, StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import '../../mol-plugin-ui/skin/light.scss';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
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

async function createBindingSiteRepresentation(plugin: PluginContext, interactions: StructureInteractions, receptors: Map<string, Structure>) {
    const contactBundles = getBindingSiteBundles(interactions, receptors);
    const update = plugin.build();

    for (const [ref, bundle] of contactBundles) {
        update.to(ref)
            .apply(StructureSelectionFromBundle, { bundle, label: 'Binding Site' })
            .apply(StructureRepresentation3D, {
                type: { name: 'ball-and-stick', params: {} },
                colorTheme: { name: 'element-symbol', params: { carbonColor: { name: 'element-symbol', params: {} } } },
            });
    }

    await update.commit();
}

function getBindingSiteBundles(interactions: StructureInteractions, receptors: Map<string, Structure>) {
    const residueIndices = new Map<string, Set<ResidueIndex>>();

    const loc = StructureElement.Location.create();
    const add = (ref: string, loci: StructureElement.Loci) => {
        if (!receptors.has(ref)) return;

        let set: Set<ResidueIndex>;
        if (residueIndices.has(ref)) {
            set = residueIndices.get(ref)!;
        } else {
            set = new Set<ResidueIndex>();
            residueIndices.set(ref, set);
        }
        StructureElement.Loci.forEachLocation(loci, l => {
            set.add(StructureProperties.residue.key(l));
        }, loc);
    };

    for (const e of interactions.elements) {
        add(e.aStructureRef!, e.a);
        add(e.bStructureRef!, e.b);
    }

    const bundles: [ref: string, bundle: StructureElement.Bundle][] = [];

    for (const [ref, indices] of Array.from(residueIndices.entries())) {
        if (indices.size === 0) continue;

        const loci = StructureQuery.loci(
            atoms({
                residueTest: e => indices.has(StructureProperties.residue.key(e.element))
            }),
            receptors.get(ref)!,
        );
        if (StructureElement.Loci.isEmpty(loci)) continue;
        bundles.push([ref, StructureElement.Bundle.fromLoci(loci)]);
    }

    return bundles;
}



async function loadComputedExample(plugin: PluginContext) {
    await plugin.clear();

    // Set up the receptor and ligand structures
    const receptorData = await plugin.builders.data.download({ url: '../../../examples/ace2.pdbqt' });
    const receptorTrajectory = await plugin.builders.structure.parseTrajectory(receptorData, 'pdbqt');
    const receptor = await plugin.builders.structure.hierarchy.applyPreset(receptorTrajectory, 'default');

    const ligandData = await plugin.builders.data.download({ url: '../../../examples/ace2-hit.mol2' });
    const ligandTrajectory = await plugin.builders.structure.parseTrajectory(ligandData, 'mol2');
    const ligand = await plugin.builders.structure.hierarchy.applyPreset(ligandTrajectory, 'default', { representationPreset: 'atomic-detail' });

    // Compute the interactions
    const update = plugin.build();

    const receptorRef = receptor?.representation.components.polymer.ref!;
    const ligandRef = ligand?.representation.components.all.ref!;

    const refs = [receptorRef, ligandRef];

    const interactionsRef = update.toRoot().apply(ComputeInteractions, {
        sources: refs.map(structureRef => ({ structureRef })),
    }, { dependsOn: refs });

    interactionsRef.apply(InteractionsShape).apply(ShapeRepresentation3D);

    await update.commit();

    console.log('Interactions', interactionsRef.selector.data?.interactions);

    // Create ball and stick representations for the binding site and focus on the ligand
    await createBindingSiteRepresentation(
        plugin,
        interactionsRef.selector.data?.interactions!,
        new Map([[receptorRef, receptor?.representation.components.polymer.data]])
    );
    PluginCommands.Camera.FocusObject(plugin, {
        targets: [{
            targetRef: ligand?.representation.representations.all.ref
        }]
    });
}

async function loadCustomExample(plugin: PluginContext) {
    await plugin.clear();

    // Set up the receptor and ligand structures
    const receptorData = await plugin.builders.data.download({ url: '../../../examples/ace2.pdbqt' });
    const receptorTrajectory = await plugin.builders.structure.parseTrajectory(receptorData, 'pdbqt');
    const receptor = await plugin.builders.structure.hierarchy.applyPreset(receptorTrajectory, 'default');

    const ligandData = await plugin.builders.data.download({ url: '../../../examples/ace2-hit.mol2' });
    const ligandTrajectory = await plugin.builders.structure.parseTrajectory(ligandData, 'mol2');
    const ligand = await plugin.builders.structure.hierarchy.applyPreset(ligandTrajectory, 'default', { representationPreset: 'atomic-detail' });

    // Compute the interactions
    const update = plugin.build();

    const receptorRef = receptor?.representation.components.polymer.ref!;
    const ligandRef = ligand?.representation.components.all.ref!;

    const refs = [receptorRef, ligandRef];

    const interactionsRef = update.toRoot().apply(CustomInteractions, {
        interactions: [
            {
                kind: 'hydrogen-bond',
                aStructureRef: receptorRef,
                a: { auth_seq_id: 353, auth_atom_id: 'N' },
                bStructureRef: ligandRef,
                b: { atom_index: 9 },
            },
            {
                kind: 'unknown',
                aStructureRef: receptorRef,
                a: { auth_seq_id: 354, auth_atom_id: 'N' },
                bStructureRef: ligandRef,
                b: [{ atom_index: 25 }, { atom_index: 26 }],
                description: 'Random test interaction'
            }
        ]
    }, { dependsOn: refs });

    interactionsRef.apply(InteractionsShape).apply(ShapeRepresentation3D);

    await update.commit();

    console.log('Interactions', interactionsRef.selector.data?.interactions);

    // Create ball and stick representations for the binding site and focus on the ligand
    await createBindingSiteRepresentation(
        plugin,
        interactionsRef.selector.data?.interactions!,
        new Map([[receptorRef, receptor?.representation.components.polymer.data]])
    );
    PluginCommands.Camera.FocusObject(plugin, {
        targets: [{
            targetRef: ligand?.representation.representations.all.ref
        }]
    });
}

const Examples = {
    computed: loadComputedExample,
    custom: loadCustomExample
};

async function init(elem: HTMLElement | string, defaultExample: keyof typeof Examples = 'computed') {
    const root = typeof elem === 'string' ? document.getElementById('viewer')! : elem;
    const plugin = await createViewer(root);
    Examples[defaultExample](plugin);
    return {
        plugin,
        examples: Object.fromEntries(Object.entries(Examples).map(([k, v]) => [k, () => v(plugin)])),
    };
}


(window as any).initInteractionsExample = init;
