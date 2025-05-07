/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRoot } from 'react-dom/client';
import { BehaviorSubject } from 'rxjs';
import { JSONCifLigandGraph, LigandGraphBondProps } from '../../extensions/json-cif/ligand-graph';
import { JSONCifDataBlock, JSONCifFile } from '../../extensions/json-cif/model';
import { ParseJSONCifFileData } from '../../extensions/json-cif/transformers';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { StructureElement, StructureProperties } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ModelFromTrajectory, StructureFromModel, TrajectoryFromMmCif } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { createPluginUI } from '../../mol-plugin-ui';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import '../../mol-plugin-ui/skin/light.scss';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { PluginSpec } from '../../mol-plugin/spec';
import { StateObjectSelector } from '../../mol-state';
import { TopologyEdits } from './edits';
import { ExampleMol } from './example-data';
import './index.html';
import { RGroupName } from './r-groups';
import { molfileToJSONCif } from './utils';

async function createViewer(root: HTMLElement) {
    const spec = DefaultPluginUISpec();
    const plugin = await createPluginUI({
        target: root,
        render: renderReact18,
        spec: {
            ...spec,
            layout: {
                initial: {
                    isExpanded: false,
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
    plugin.managers.interactivity.setProps({ granularity: 'element' });
    plugin.selectionMode = true;

    return plugin;
}

async function loadMolfile(model: EditorModel, molfile: string) {
    const { plugin } = model;

    await plugin.clear();

    const file = await molfileToJSONCif(molfile);

    const update = plugin.build();

    const data = update.toRoot()
        .apply(ParseJSONCifFileData, { data: file.data });

    data
        .apply(TrajectoryFromMmCif)
        .apply(ModelFromTrajectory)
        .apply(StructureFromModel, { type: { name: 'model', params: { } } })
        .apply(StructureRepresentation3D, {
            type: { name: 'ball-and-stick', params: { } },
            colorTheme: {
                name: 'element-symbol',
                params: { carbonColor: { name: 'element-symbol', params: {} } }
            }
        });

    await update.commit();

    model.dataSelector = data.selector;
}

class EditorModel {
    dataSelector: StateObjectSelector | undefined = undefined;
    history: JSONCifFile[] = [];

    state = {
        element: new BehaviorSubject<string>('C'),
    };

    get data() {
        return this.dataSelector?.cell?.transform?.params?.data as JSONCifFile | undefined;
    }

    createGraph() {
        return new JSONCifLigandGraph(this.data?.dataBlocks[0]!);
    }

    async update(data: JSONCifDataBlock) {
        if (!this.data) return;

        const updated: JSONCifFile = {
            ...this.data!,
            dataBlocks: [data],
        };

        this.history.push(this.data!);
        const update = this.plugin.build();
        update.to(this.dataSelector!).update({ data: updated });
        await update.commit();
    }

    undo = async () => {
        if (!this.dataSelector) return;
        if (this.history.length === 0) return;

        const data = this.history.pop()!;
        const update = this.plugin.build();
        update.to(this.dataSelector).update({ data });
        await update.commit();
    };

    private getEditableStructures() {
        if (!this.dataSelector?.isOk) return new Set();

        const structures = this.plugin.state.data.selectQ(q => q
            .byRef(this.dataSelector?.ref!)
            .subtree()
            .filter(c => PluginStateObject.Molecule.Structure.is(c.obj))
        );
        return new Set(structures.map(s => s.obj?.data));
    }

    private getSelectedAtomIds() {
        if (!this.data) return [];

        const structures = this.getEditableStructures();
        if (structures.size === 0) return [];

        const { selection } = this.plugin.managers.structure;
        const ids: number[] = [];
        selection.entries.forEach(e => {
            if (!structures.has(e.selection.structure)) return;
            StructureElement.Loci.forEachLocation(e.selection, (l) => {
                ids.push(StructureProperties.atom.id(l));
            });
        });
        return ids;
    }

    async editGraphTopology<Args extends any[], T>(fn: (graph: JSONCifLigandGraph, ...args: Args) => Promise<T>, ...args: Args) {
        try {
            const graph = this.createGraph();
            const result = await fn(graph, ...args);
            const data = graph.getData().block;
            await this.update(data);
            return result;
        } catch (e) {
            console.error('Failed to edit graph');
            console.error(e);
        }
    }

    setElement = async () => {
        const symbol = this.state.element.value.trim();
        if (!symbol) return;

        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        await this.editGraphTopology(TopologyEdits.setElement, ids, symbol);
    };

    addElement = async () => {
        const symbol = this.state.element.value.trim();
        if (!symbol) return;

        const ids = this.getSelectedAtomIds();
        if (ids.length !== 1) return;

        await this.editGraphTopology(TopologyEdits.addElement, ids[0], symbol);
    };

    removeAtoms = async () => {
        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        await this.editGraphTopology(TopologyEdits.removeAtoms, ids);
    };

    removeBonds = async () => {
        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        await this.editGraphTopology(TopologyEdits.removeBonds, ids);
    };

    updateBonds = async (props: LigandGraphBondProps) => {
        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        await this.editGraphTopology(TopologyEdits.updateBonds, ids, props);
    };

    attachRgroup = async (name: RGroupName) => {
        const ids = this.getSelectedAtomIds();
        if (ids.length !== 1) return;

        await this.editGraphTopology(TopologyEdits.attachRgroup, ids[0], name);
    };

    constructor(public plugin: PluginContext) { }
}

function ControlsUI({ model }: { model: EditorModel }) {
    return <div>
        <EditElementSymbolUI model={model} />
    </div>;
}

function EditElementSymbolUI({ model }: { model: EditorModel }) {
    return <div style={{ display: 'flex', flexDirection: 'column', gap: '5px' }}>
        <div style={{ display: 'flex', gap: '5px' }}>
            <b>Atoms:</b>
            <button onClick={model.removeAtoms}>Remove</button>
            <div>
                <ElementEditUI model={model} />
                <button onClick={model.setElement}>Set</button>
                <button onClick={model.addElement}>Add</button>
            </div>
        </div>
        <div style={{ display: 'flex', gap: '5px' }}>
            <b>Bonds:</b>
            <button onClick={model.removeBonds}>Remove</button>
            <button onClick={() => model.updateBonds({ value_order: 'sing', type_id: 'covale' })}>-</button>
            <button onClick={() => model.updateBonds({ value_order: 'doub', type_id: 'covale' })}>=</button>
            <button onClick={() => model.updateBonds({ value_order: 'trip', type_id: 'covale' })}>≡</button>
        </div>
        <div style={{ display: 'flex', gap: '5px' }}>
            <b>R-groups:</b>
            <button onClick={() => model.attachRgroup('CH3')}>CH3</button>
        </div>
        <div>
            <button onClick={model.undo}>Undo</button>
        </div>
    </div>;
}

function ElementEditUI({ model }: { model: EditorModel }) {
    const element = useBehavior(model.state.element);
    return <input type="text" value={element} style={{ width: 50 }} onChange={e => model.state.element.next(e.target.value)} />;
}

async function init(viewer: HTMLElement | string, controls: HTMLElement | string) {
    const root = typeof viewer === 'string' ? document.getElementById('viewer')! : viewer;
    const plugin = await createViewer(root);

    const model = new EditorModel(plugin);

    createRoot(
        typeof controls === 'string' ? document.getElementById('controls')! : controls
    ).render(<ControlsUI model={model} />);

    loadMolfile(model, ExampleMol);
    return model;
}

(window as any).initLigandEditorExample = init;
