/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRoot } from 'react-dom/client';
import { BehaviorSubject } from 'rxjs';
import { InteractionElementSchema, InteractionKind, StructureInteractionElement, StructureInteractions } from '../../extensions/interactions/model';
import { ComputeContacts, CustomInteractions, InteractionsShape } from '../../extensions/interactions/transforms';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { ResidueIndex, Structure, StructureElement, StructureProperties, StructureQuery } from '../../mol-model/structure';
import { atoms } from '../../mol-model/structure/query/queries/generators';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { MultiStructureSelectionFromBundle, StructureSelectionFromBundle } from '../../mol-plugin-state/transforms/model';
import { ShapeRepresentation3D, StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { createPluginUI } from '../../mol-plugin-ui';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { PluginSpec } from '../../mol-plugin/spec';

import '../../mol-plugin-ui/skin/light.scss';
import './index.html';
import { Task } from '../../mol-task';
import { computeContacts } from '../../extensions/interactions/compute';

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
            ]
        }
    });

    return plugin;
}

async function createBindingSiteRepresentation(plugin: PluginContext, interactions: StructureInteractions[], receptors: Map<string, Structure>) {
    const contactBundles = getBindingSiteBundles(interactions.flatMap(e => e.elements), receptors);
    const update = plugin.build();

    for (const [ref, bundle] of contactBundles) {
        update.to(ref)
            .apply(StructureSelectionFromBundle, { bundle, label: 'Binding Site' })
            .apply(StructureRepresentation3D, {
                type: { name: 'ball-and-stick', params: { sizeFactor: 0.2 } },
                colorTheme: { name: 'element-symbol', params: { carbonColor: { name: 'element-symbol', params: {} } } },
            });
    }

    await update.commit();
}

function getBindingSiteBundles(interactions: StructureInteractionElement[], receptors: Map<string, Structure>) {
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

    for (const e of interactions) {
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

async function loadComputedExample(
    plugin: PluginContext,
    { receptorUrl, ligandUrl }: { receptorUrl: [url: string, format: BuiltInTrajectoryFormat], ligandUrl: [url: string, format: BuiltInTrajectoryFormat] },
    options: { receptor_label_asym_id: string | undefined, analyzeTrajectory?: boolean }
) {
    await plugin.clear();

    // Set up the receptor and ligand structures
    const receptorData = await plugin.builders.data.download({ url: receptorUrl[0] });
    const receptorTrajectory = await plugin.builders.structure.parseTrajectory(receptorData, receptorUrl[1]);
    const receptor = await plugin.builders.structure.hierarchy.applyPreset(receptorTrajectory, 'default', { representationPreset: 'polymer-cartoon' });

    const ligandData = await plugin.builders.data.download({ url: ligandUrl[0] });
    const ligandTrajectory = await plugin.builders.structure.parseTrajectory(ligandData, ligandUrl[1]);
    const ligand = await plugin.builders.structure.hierarchy.applyPreset(ligandTrajectory, 'default', { representationPreset: 'atomic-detail' });

    // Compute the interactions
    const update = plugin.build();

    const receptorRef = receptor?.structure.ref!;
    const ligandRef = ligand?.structure.ref!;

    const refs = [receptorRef, ligandRef];
    const interactionsRef = update.toRoot()
        .apply(MultiStructureSelectionFromBundle, {
            selections: [
                { key: 'a', ref: receptorRef, bundle: StructureElement.Schema.toBundle(receptor?.structure.data!, { label_asym_id: options.receptor_label_asym_id }) },
                { key: 'b', ref: ligandRef, bundle: StructureElement.Schema.toBundle(ligand?.structure.data!, { }) },
            ],
            isTransitive: true,
            label: 'Label'
        }, { dependsOn: refs })
        .apply(ComputeContacts);

    interactionsRef.apply(InteractionsShape).apply(ShapeRepresentation3D);

    await update.commit();

    if (!options.analyzeTrajectory) {
        console.log('Interactions', interactionsRef.selector.data?.interactions);

        // Create ball and stick representations for the binding site and focus on the ligand
        await createBindingSiteRepresentation(
            plugin,
            [interactionsRef.selector.data?.interactions!],
            new Map([[receptorRef, receptor?.structure.data!]])
        );
    } else {
        const trajectoryInteractions: StructureInteractions[] = [];
        const receptorLoci = StructureElement.Schema.toLoci(receptor?.structure.data!, { label_asym_id: options.receptor_label_asym_id });
        for (let fI = 0; fI < ligandTrajectory.data!.frameCount; fI++) {
            const model = await Task.resolveInContext(ligandTrajectory.data!.getFrameAtIndex(fI));
            const structure = Structure.ofModel(model);
            const currentInteractions = await plugin.runTask(Task.create('Compute Contacts', ctx => {
                return computeContacts(ctx, [
                    { structureRef: receptorRef, loci: receptorLoci },
                    { structureRef: ligandRef, loci: Structure.toStructureElementLoci(structure) },
                ]);
            }));
            trajectoryInteractions.push(currentInteractions);
        }

        console.log('Interactions', trajectoryInteractions);

        await createBindingSiteRepresentation(
            plugin,
            trajectoryInteractions,
            new Map([[receptorRef, receptor?.structure.data!]])
        );
    }

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
            }
        ]
    }, { dependsOn: refs });

    interactionsRef.apply(InteractionsShape).apply(ShapeRepresentation3D);

    await update.commit();

    console.log('Interactions', interactionsRef.selector.data?.interactions);

    // Create ball and stick representations for the binding site and focus on the ligand
    await createBindingSiteRepresentation(
        plugin,
        [interactionsRef.selector.data?.interactions!],
        new Map([[receptorRef, receptor?.representation.components.polymer.data]])
    );
    PluginCommands.Camera.FocusObject(plugin, {
        targets: [{
            targetRef: ligand?.representation.representations.all.ref
        }]
    });
}

async function loadTestAllExample(plugin: PluginContext) {
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

    const basic = (kind: InteractionKind, atom_index: number | number[], description?: string): InteractionElementSchema => {
        return {
            kind,
            aStructureRef: receptorRef,
            a: { auth_seq_id: 354, auth_atom_id: 'N' },
            bStructureRef: ligandRef,
            b: Array.isArray(atom_index) ? { items: { atom_index } } : { atom_index },
            description,
        };
    };

    const covalent = (degree: number, atom_index: number): InteractionElementSchema => {
        return {
            kind: 'covalent',
            degree: degree === -1 ? 'aromatic' : Math.abs(degree) as 1 | 2 | 3 | 4,
            aStructureRef: receptorRef,
            a: { auth_seq_id: 354, auth_atom_id: 'N' },
            bStructureRef: ligandRef,
            b: { atom_index }
        };
    };

    const interactionsRef = update.toRoot().apply(CustomInteractions, {
        interactions: [
            basic('unknown', 1),
            basic('ionic', 2),
            basic('pi-stacking', 3),
            basic('cation-pi', 4),
            basic('halogen-bond', 5),
            basic('hydrogen-bond', 6),
            basic('weak-hydrogen-bond', 7),
            basic('hydrophobic', 8),
            basic('metal-coordination', 9),
            covalent(1, 10),
            covalent(2, 11),
            covalent(3, 12),
            covalent(-1, 13), // aromatic
            basic('unknown', [0, 1, 2, 3, 13, 14], 'Testing centroid for atom set'),
        ]
    }, { dependsOn: refs });

    interactionsRef.apply(InteractionsShape).apply(ShapeRepresentation3D);

    await update.commit();

    console.log('Interactions', interactionsRef.selector.data?.interactions);

    // Create ball and stick representations for the binding site and focus on the ligand
    await createBindingSiteRepresentation(
        plugin,
        [interactionsRef.selector.data?.interactions!],
        new Map([[receptorRef, receptor?.representation.components.polymer.data]])
    );
    PluginCommands.Camera.FocusObject(plugin, {
        targets: [{
            targetRef: ligand?.representation.representations.all.ref
        }]
    });
}

const Examples = {
    'Computed (1iep)': (plugin: PluginContext) => loadComputedExample(plugin, {
        receptorUrl: ['https://files.rcsb.org/download/1IEP.cif', 'mmcif'],
        ligandUrl: ['https://models.rcsb.org/v1/1iep/atoms?label_asym_id=G&copy_all_categories=false', 'mmcif']
    }, { receptor_label_asym_id: 'A' }),
    'Computed (ACE2)': (plugin: PluginContext) => loadComputedExample(plugin, {
        receptorUrl: ['../../../examples/ace2.pdbqt', 'pdbqt'],
        ligandUrl: ['../../../examples/ace2-hit.mol2', 'mol2']
    }, { receptor_label_asym_id: 'B' }),
    'Computed (multiple)': (plugin: PluginContext) => loadComputedExample(plugin, {
        receptorUrl: ['../../../examples/docking/receptor_1.pdb', 'pdb'],
        ligandUrl: ['../../../examples/docking/ligands_1.sdf', 'sdf']
    }, { receptor_label_asym_id: undefined, analyzeTrajectory: true }),
    'Custom': loadCustomExample,
    'Synthetic': loadTestAllExample
};

function SelectExampleUI({ state, load }: {
    state: BehaviorSubject<{ name?: keyof typeof Examples, isLoading?: boolean }>,
    load: (name: keyof typeof Examples) => void
}) {
    const current = useBehavior(state);
    return <div>
        Select Example:{' '}
        <select value={current.name} onChange={e => load(e.target.value as any)} disabled={current.isLoading}>
            {Object.keys(Examples).map(k => <option key={k} value={k}>{k}</option>)}
        </select>
    </div>;
}

async function init(viewer: HTMLElement | string, controls: HTMLElement | string, defaultExample: keyof typeof Examples = 'Computed (1iep)') {
    const root = typeof viewer === 'string' ? document.getElementById(viewer)! : viewer;
    const plugin = await createViewer(root);

    const state = new BehaviorSubject<{ name?: keyof typeof Examples, isLoading?: boolean }>({});
    const loadExample = async (name: keyof typeof Examples) => {
        state.next({ name, isLoading: true });
        try {
            await Examples[name](plugin);
            state.next({ name });
        } catch (e) {
            console.error(e);
            state.next({});
        }
    };

    createRoot(
        typeof controls === 'string' ? document.getElementById(controls)! : controls
    ).render(<SelectExampleUI state={state} load={loadExample} />);

    loadExample(defaultExample);
    return { plugin, loadExample };
}


(window as any).initInteractionsExample = init;
