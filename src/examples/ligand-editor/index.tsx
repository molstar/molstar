/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRoot } from 'react-dom/client';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { ModelFromTrajectory, StructureFromModel, TrajectoryFromMmCif } from '../../mol-plugin-state/transforms/model';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { PluginSpec } from '../../mol-plugin/spec';
import { ParseJSONCifFileData } from '../../extensions/json-cif/transformers';
import '../../mol-plugin-ui/skin/light.scss';
import './index.html';
import { getJSONCifFile } from './load';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { StateObjectSelector } from '../../mol-state';
import { JSONCifDataBlock, JSONCifFile } from '../../extensions/json-cif/model';
import { StructureElement, StructureProperties } from '../../mol-model/structure';
import { JSONCifLigandGraph, LigandGraphBondProps } from '../../extensions/json-cif/ligand-graph';
import { BehaviorSubject } from 'rxjs';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';

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

async function loadMolfile(model: EditorModel, molfile?: string) {
    const { plugin } = model;

    await plugin.clear();

    const file = await getJSONCifFile(molfile);

    const update = plugin.build();

    const data = update.toRoot()
        .apply(ParseJSONCifFileData, { data: file });

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

    private getSelectedAtomIds() {
        const selection = this.plugin.managers.structure.selection;
        const ids: number[] = [];
        selection.entries.forEach(e => {
            StructureElement.Loci.forEachLocation(e.selection, (l) => {
                ids.push(StructureProperties.atom.id(l));
            });
        });
        return ids;
    }

    setElement = async () => {
        const symbol = this.state.element.value.trim();
        if (!symbol) return;

        const { data } = this;
        if (!data) return;

        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        const graph = this.createGraph();
        for (const id of ids) {
            graph.modifyAtom(id, { type_symbol: symbol });
        }

        await this.update(graph.getData().block);
    };

    addElement = async () => {
        const symbol = this.state.element.value.trim();
        if (!symbol) return;

        const { data } = this;
        if (!data) return;

        const ids = this.getSelectedAtomIds();
        if (ids.length !== 1) return;

        const graph = this.createGraph();
        for (const id of ids) {
            graph.attachAtom(id, { type_symbol: symbol });
        }

        await this.update(graph.getData().block);
    };

    removeAtoms = async () => {
        const { data } = this;
        if (!data) return;

        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        const graph = this.createGraph();
        for (const id of ids) {
            graph.removeAtom(id);
        }

        await this.update(graph.getData().block);
    };

    removeBonds = async () => {
        const { data } = this;
        if (!data) return;

        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        const graph = this.createGraph();

        for (let i = 0; i < ids.length; ++i) {
            for (let j = i + 1; j < ids.length; ++j) {
                graph.removeBond(ids[i], ids[j]);
            }
        }

        await this.update(graph.getData().block);
    };

    updateBonds = async (props: LigandGraphBondProps) => {
        const { data } = this;
        if (!data) return;

        const ids = this.getSelectedAtomIds();
        if (!ids.length) return;

        const graph = this.createGraph();

        // TODO: iterate on the all-pairs behavior
        // e.g. only add bonds if there is no path connecting them,
        // or by a distance threshold, ...
        for (let i = 0; i < ids.length; ++i) {
            for (let j = i + 1; j < ids.length; ++j) {
                graph.addOrUpdateBond(ids[i], ids[j], props);
            }
        }

        await this.update(graph.getData().block);
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
            <div>TODO</div>
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

    loadMolfile(model);
    return model;
}

(window as any).initLigandEditorExample = init;
