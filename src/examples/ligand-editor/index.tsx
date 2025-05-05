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

import { useState } from 'react';
import { ParseJSONCifFileData } from '../../extensions/json-cif/transformers';
import '../../mol-plugin-ui/skin/light.scss';
import './index.html';
import { getJSONCifFile } from './load';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { StateObjectSelector } from '../../mol-state';
import { JSONCifFile } from '../../extensions/json-cif/model';
import { StructureElement, StructureProperties } from '../../mol-model/structure';
import { produce } from 'immer';

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
        .apply(StructureFromModel, { type: { name: 'model', params: {} } })
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

    get data() {
        return this.dataSelector?.cell?.transform?.params?.data as JSONCifFile | undefined;
    }

    async updateData(data: JSONCifFile) {
        if (!this.dataSelector) return;

        const update = this.plugin.build();
        update.to(this.dataSelector).update({ data });
        await update.commit();
    }

    async setElementSymbol(symbol: string) {
        const { data } = this;
        if (!data) return;

        const selection = this.plugin.managers.structure.selection;
        let fst: StructureElement.Location | undefined;
        selection.entries.forEach(e => {
            if (fst) return;
            fst = StructureElement.Loci.getFirstLocation(e.selection);
        });

        if (!fst) return;

        const index = StructureProperties.atom.sourceIndex(fst);

        const updatedData = produce(data, draft => {
            draft.dataBlocks[0].categories['atom_site'].rows[index]['type_symbol'] = symbol;
        });

        await this.updateData(updatedData);
    }

    constructor(public plugin: PluginContext) { }
}

function ControlsUI({ model }: { model: EditorModel }) {
    return <div>
        <EditElementSymbolUI model={model} />
    </div>;
}

function EditElementSymbolUI({ model }: { model: EditorModel }) {
    const [symbol, setSymbol] = useState('C');

    return <div>
        <input type="text" value={symbol} onChange={e => setSymbol(e.target.value)} />
        <button onClick={() => model.setElementSymbol(symbol)}>Set Element Symbol</button>
    </div>;
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
