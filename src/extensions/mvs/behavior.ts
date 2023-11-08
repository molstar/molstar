/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomModelProperty } from '../../mol-model-props/common/custom-model-property';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { LociLabelProvider } from '../../mol-plugin-state/manager/loci-label';
import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { StructureRepresentationProvider } from '../../mol-repr/structure/representation';
import { ColorTheme } from '../../mol-theme/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { AnnotationColorThemeProvider } from './additions/annotation-color-theme';
import { AnnotationLabelRepresentationProvider } from './additions/annotation-label/representation';
import { AnnotationsProvider } from './additions/annotation-prop';
import { AnnotationTooltipsLabelProvider, AnnotationTooltipsProvider } from './additions/annotation-tooltips-prop';
import { CustomLabelRepresentationProvider } from './additions/custom-label/representation';
import { CustomTooltipsLabelProvider, CustomTooltipsProvider } from './additions/custom-tooltips-prop';
import { makeMultilayerColorThemeProvider } from './additions/multilayer-color-theme';


/** Collection of things that can be register/unregistered in a plugin */
interface Registrables {
    customModelProperties?: CustomModelProperty.Provider<any, any>[],
    customStructureProperties?: CustomStructureProperty.Provider<any, any>[],
    representations?: StructureRepresentationProvider<any>[],
    colorThemes?: ColorTheme.Provider[],
    lociLabels?: LociLabelProvider[],
}


/** Registers everything needed for loading MolViewSpec files */
export const MolViewSpec = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'molviewspec',
    category: 'misc',
    display: {
        name: 'MolViewSpec',
        description: 'MolViewSpec extension'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private readonly registrables: Registrables = {
            customModelProperties: [
                AnnotationsProvider,
            ],
            customStructureProperties: [
                CustomTooltipsProvider,
                AnnotationTooltipsProvider,
            ],
            representations: [
                CustomLabelRepresentationProvider,
                AnnotationLabelRepresentationProvider,
            ],
            colorThemes: [
                AnnotationColorThemeProvider,
                makeMultilayerColorThemeProvider(this.ctx.representation.structure.themes.colorThemeRegistry),
            ],
            lociLabels: [
                CustomTooltipsLabelProvider,
                AnnotationTooltipsLabelProvider,
            ],
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
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
    })
});
