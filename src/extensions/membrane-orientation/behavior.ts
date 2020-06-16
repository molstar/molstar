/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '../../mol-plugin-state/builder/structure/representation-preset';
import { MembraneOrientationProvider, MembraneOrientation } from './membrane-orientation';
import { StateObjectRef } from '../../mol-state';
import { Task } from '../../mol-task';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { MembraneOrientationRepresentationProvider } from './representation';
import { HydrophobicityColorThemeProvider } from '../../mol-theme/color/hydrophobicity';

export const MembraneOrientationData = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'membrane-orientation-prop',
    category: 'custom-props',
    display: {
        name: 'Membrane Orientation',
        description: 'Initialize orientation of membrane layers. Data calculated with ANVIL algorithm.' // TODO add ' or obtained via RCSB PDB'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = MembraneOrientationProvider

        register(): void {
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

            this.ctx.representation.structure.themes.colorThemeRegistry.add(HydrophobicityColorThemeProvider);

            this.ctx.representation.structure.registry.add(MembraneOrientationRepresentationProvider);
            this.ctx.query.structure.registry.add(MembraneOrientation.isTransmembrane);

            this.ctx.builders.structure.representation.registerPreset(MembraneOrientationPreset);
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);

            this.ctx.representation.structure.themes.colorThemeRegistry.remove(HydrophobicityColorThemeProvider);

            this.ctx.representation.structure.registry.remove(MembraneOrientationRepresentationProvider);
            this.ctx.query.structure.registry.remove(MembraneOrientation.isTransmembrane);

            this.ctx.builders.structure.representation.unregisterPreset(MembraneOrientationPreset);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});

export const MembraneOrientationPreset = StructureRepresentationPresetProvider({
    id: 'preset-membrane-orientation',
    display: {
        name: 'Membrane Orientation', group: 'Annotation',
        description: 'Initialize orientation of membrane layers. Data calculated with ANVIL algorithm.' // TODO add ' or obtained via RCSB PDB'
    },
    isApplicable(a) {
        return MembraneOrientationProvider.isApplicable(a.data);
    },
    params: () => StructureRepresentationPresetProvider.CommonParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const structure  = structureCell?.obj?.data;
        if (!structureCell || !structure) return {};

        await plugin.runTask(Task.create('Membrane Orientation', async runtime => {
            await MembraneOrientationProvider.attach({ runtime, assetManager: plugin.managers.asset }, structure);
        }));

        const colorTheme = HydrophobicityColorThemeProvider.name as any;
        return await PresetStructureRepresentations.auto.apply(ref, { ...params, theme: { globalName: colorTheme, focus: { name: colorTheme } } }, plugin);
    }
});
