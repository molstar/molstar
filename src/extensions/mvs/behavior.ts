/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomModelProperty } from '../../mol-model-props/common/custom-model-property';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { DataFormatProvider } from '../../mol-plugin-state/formats/provider';
import { PluginDragAndDropHandler } from '../../mol-plugin-state/manager/drag-and-drop';
import { LociLabelProvider } from '../../mol-plugin-state/manager/loci-label';
import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { StructureRepresentationProvider } from '../../mol-repr/structure/representation';
import { StateAction, StateObjectCell, StateTree } from '../../mol-state';
import { Task } from '../../mol-task';
import { ColorTheme } from '../../mol-theme/color';
import { fileToDataUri } from '../../mol-util/file';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { MVSAnnotationColorThemeProvider } from './components/annotation-color-theme';
import { MVSAnnotationLabelRepresentationProvider } from './components/annotation-label/representation';
import { MVSAnnotationsProvider } from './components/annotation-prop';
import { MVSAnnotationTooltipsLabelProvider, MVSAnnotationTooltipsProvider } from './components/annotation-tooltips-prop';
import { CustomLabelRepresentationProvider } from './components/custom-label/representation';
import { CustomTooltipsLabelProvider, CustomTooltipsProvider } from './components/custom-tooltips-prop';
import { LoadMvsData, MVSJFormatProvider, MVSXFormatProvider, loadMVSX } from './components/formats';
import { IsMVSModelProvider } from './components/is-mvs-model-prop';
import { makeMultilayerColorThemeProvider } from './components/multilayer-color-theme';
import { loadMVS } from './load';
import { MVSData } from './mvs-data';


/** Collection of things that can be register/unregistered in a plugin */
interface Registrables {
    customModelProperties?: CustomModelProperty.Provider<any, any>[],
    customStructureProperties?: CustomStructureProperty.Provider<any, any>[],
    representations?: StructureRepresentationProvider<any>[],
    colorThemes?: ColorTheme.Provider[],
    lociLabels?: LociLabelProvider[],
    dragAndDropHandlers?: DragAndDropHandler[],
    dataFormats?: { name: string, provider: DataFormatProvider }[],
    actions?: StateAction[],
}


/** Registers everything needed for loading MolViewSpec files */
export const MolViewSpec = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'molviewspec',
    category: 'misc',
    display: {
        name: 'MolViewSpec',
        description: 'MolViewSpec extension',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private readonly registrables: Registrables = {
            customModelProperties: [
                IsMVSModelProvider,
                MVSAnnotationsProvider,
            ],
            customStructureProperties: [
                CustomTooltipsProvider,
                MVSAnnotationTooltipsProvider,
            ],
            representations: [
                CustomLabelRepresentationProvider,
                MVSAnnotationLabelRepresentationProvider,
            ],
            colorThemes: [
                MVSAnnotationColorThemeProvider,
                makeMultilayerColorThemeProvider(this.ctx.representation.structure.themes.colorThemeRegistry),
            ],
            lociLabels: [
                CustomTooltipsLabelProvider,
                MVSAnnotationTooltipsLabelProvider,
            ],
            dragAndDropHandlers: [
                MVSDragAndDropHandler,
            ],
            dataFormats: [
                { name: 'MVSJ', provider: MVSJFormatProvider },
                { name: 'MVSX', provider: MVSXFormatProvider },
            ],
            actions: [
                LoadMvsData,
            ]
        };

        register(): void {
            for (const prop of this.registrables.customModelProperties ?? []) {
                this.ctx.customModelProperties.register(prop, this.params.autoAttach);
            }
            for (const prop of this.registrables.customStructureProperties ?? []) {
                this.ctx.customStructureProperties.register(prop, this.params.autoAttach);
            }
            for (const repr of this.registrables.representations ?? []) {
                this.ctx.representation.structure.registry.add(repr);
            }
            for (const theme of this.registrables.colorThemes ?? []) {
                this.ctx.representation.structure.themes.colorThemeRegistry.add(theme);
            }
            for (const provider of this.registrables.lociLabels ?? []) {
                this.ctx.managers.lociLabels.addProvider(provider);
            }
            for (const handler of this.registrables.dragAndDropHandlers ?? []) {
                this.ctx.managers.dragAndDrop.addHandler(handler.name, handler.handle);
            }
            for (const format of this.registrables.dataFormats ?? []) {
                this.ctx.dataFormats.add(format.name, format.provider);
            }
            for (const action of this.registrables.actions ?? []) {
                this.ctx.state.data.actions.add(action);
            }

            this.ctx.managers.markdownExtensions.registerRefResolver('mvs', (plugin, refs) => {
                const mvsRefs = new Set(refs.map(ref => `mvs-ref:${ref}`));
                return StateTree.doPreOrder(
                    plugin.state.data.tree,
                    plugin.state.data.tree.root,
                    { mvsRefs, plugin, cells: [] as StateObjectCell[] },
                    (n, _, s) => {
                    if (!n.tags) return;
                    for (const tag of n.tags) {
                        if (!s.mvsRefs.has(tag)) continue;
                        const cell = s.plugin.state.data.cells.get(n.ref);
                        if (cell) {
                            s.cells.push(cell);
                            break;
                        }
                    }
                }).cells;
            });

            this.ctx.managers.markdownExtensions.registerUriResolver('mvs', (plugin, uri) => {
                const { assets } = plugin.managers.asset;
                const asset = assets.find(a => a.file.name === uri);
                if (!asset) {
                    return undefined;
                }
                try {
                    return fileToDataUri(asset.file);
                } catch (e) {
                    console.error(`MVS: Failed to convert asset file to data URI for '${uri}'`, e);
                    return undefined;
                }
            });
        }
        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            for (const prop of this.registrables.customModelProperties ?? []) {
                this.ctx.customModelProperties.setDefaultAutoAttach(prop.descriptor.name, this.params.autoAttach);
            }
            for (const prop of this.registrables.customStructureProperties ?? []) {
                this.ctx.customStructureProperties.setDefaultAutoAttach(prop.descriptor.name, this.params.autoAttach);
            }
            return updated;
        }
        unregister() {
            for (const prop of this.registrables.customModelProperties ?? []) {
                this.ctx.customModelProperties.unregister(prop.descriptor.name);
            }
            for (const prop of this.registrables.customStructureProperties ?? []) {
                this.ctx.customStructureProperties.unregister(prop.descriptor.name);
            }
            for (const repr of this.registrables.representations ?? []) {
                this.ctx.representation.structure.registry.remove(repr);
            }
            for (const theme of this.registrables.colorThemes ?? []) {
                this.ctx.representation.structure.themes.colorThemeRegistry.remove(theme);
            }
            for (const labelProvider of this.registrables.lociLabels ?? []) {
                this.ctx.managers.lociLabels.removeProvider(labelProvider);
            }
            for (const handler of this.registrables.dragAndDropHandlers ?? []) {
                this.ctx.managers.dragAndDrop.removeHandler(handler.name);
            }
            for (const format of this.registrables.dataFormats ?? []) {
                this.ctx.dataFormats.remove(format.name);
            }
            for (const action of this.registrables.actions ?? []) {
                this.ctx.state.data.actions.remove(action);
            }
            this.ctx.managers.markdownExtensions.removeRefResolver('mvs');
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
    })
});


/** Registrable method for handling dragged-and-dropped files */
interface DragAndDropHandler {
    name: string,
    handle: PluginDragAndDropHandler,
}

/** DragAndDropHandler handler for `.mvsj` and `.mvsx` files */
const MVSDragAndDropHandler: DragAndDropHandler = {
    name: 'mvs-mvsj-mvsx',
    /** Load .mvsj and .mvsx files. Delete previous plugin state before loading.
     * If multiple files are provided, merge their MVS data into one state.
     * Return `true` if at least one file has been loaded. */
    async handle(files: File[], plugin: PluginContext): Promise<boolean> {
        let applied = false;
        for (const file of files) {
            if (file.name.toLowerCase().endsWith('.mvsj')) {
                const task = Task.create('Load MVSJ file', async ctx => {
                    const data = await file.text();
                    const mvsData = MVSData.fromMVSJ(data);
                    await loadMVS(plugin, mvsData, { sanityChecks: true, appendSnapshots: applied, sourceUrl: undefined });
                });
                await plugin.runTask(task);
                applied = true;
            }
            if (file.name.toLowerCase().endsWith('.mvsx')) {
                const task = Task.create('Load MVSX file', async ctx => {
                    const buffer = await file.arrayBuffer();
                    const array = new Uint8Array(buffer);
                    const parsed = await loadMVSX(plugin, ctx, array);
                    await loadMVS(plugin, parsed.mvsData, { sanityChecks: true, appendSnapshots: applied, sourceUrl: parsed.sourceUrl });
                });
                await plugin.runTask(task);
                applied = true;
            }
        }
        return applied;
    },
};
