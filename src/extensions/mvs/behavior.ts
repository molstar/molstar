/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { StateAction } from '../../mol-state';
import { ColorTheme } from '../../mol-theme/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { MVSAnnotationColorThemeProvider } from './components/annotation-color-theme';
import { MVSAnnotationLabelRepresentationProvider } from './components/annotation-label/representation';
import { MVSAnnotationsProvider } from './components/annotation-prop';
import { MVSAnnotationTooltipsLabelProvider, MVSAnnotationTooltipsProvider } from './components/annotation-tooltips-prop';
import { CustomLabelRepresentationProvider } from './components/custom-label/representation';
import { CustomTooltipsLabelProvider, CustomTooltipsProvider } from './components/custom-tooltips-prop';
import { LoadMvsData, MVSJFormatProvider } from './components/formats';
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

/** DragAndDropHandler handler for `.mvsj` files */
const MVSDragAndDropHandler: DragAndDropHandler = {
    name: 'mvs-mvsj',
    /** Load .mvsj files. Delete previous plugin state before loading.
     * If multiple files are provided, merge their MVS data into one state. */
    async handle(files: File[], plugin: PluginContext): Promise<boolean> {
        let applied = false;
        for (const file of files) {
            if (file.name.toLowerCase().endsWith('.mvsj')) {
                const data = await file.text();
                const mvsData = MVSData.fromMVSJ(data);
                await loadMVS(plugin, mvsData, { sanityChecks: true, replaceExisting: !applied });
                applied = true;
            }
        }
        return applied;
    },
};
