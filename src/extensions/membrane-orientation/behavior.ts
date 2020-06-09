/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '../../mol-plugin-state/builder/structure/representation-preset';
import { MembraneOrientationProvider } from '../../mol-model-props/computed/membrane-orientation';
import { StateObjectRef } from '../../mol-state';
import { Task } from '../../mol-task';
import { AccessibleSurfaceAreaColorThemeProvider } from '../../mol-model-props/computed/themes/accessible-surface-area';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { StructureSelectionQuery, StructureSelectionCategory } from '../../mol-plugin-state/helpers/structure-selection-query';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { DefaultQueryRuntimeTable, QuerySymbolRuntime } from '../../mol-script/runtime/query/base';
import { MembraneOrientationRepresentationProvider } from './representation';
import { CustomPropSymbol } from '../../mol-script/language/symbol';
import Type from '../../mol-script/language/type';
import { isInMembranePlane } from '../../mol-model-props/computed/membrane-orientation/ANVIL';
import { StructureProperties } from '../../mol-model/structure';
import { Vec3 } from '../../mol-math/linear-algebra';

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
            DefaultQueryRuntimeTable.addCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

            this.ctx.representation.structure.themes.colorThemeRegistry.add(AccessibleSurfaceAreaColorThemeProvider);

            this.ctx.representation.structure.registry.add(MembraneOrientationRepresentationProvider);
            this.ctx.query.structure.registry.add(isTransmembrane);

            this.ctx.builders.structure.representation.registerPreset(MembraneOrientationPreset);
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            DefaultQueryRuntimeTable.removeCustomProp(this.provider.descriptor);
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);

            this.ctx.representation.structure.themes.colorThemeRegistry.remove(AccessibleSurfaceAreaColorThemeProvider);

            this.ctx.representation.structure.registry.remove(MembraneOrientationRepresentationProvider);
            this.ctx.query.structure.registry.remove(isTransmembrane);

            this.ctx.builders.structure.representation.unregisterPreset(MembraneOrientationPreset);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});

const pos = Vec3();
const isTransmembraneSymbol = QuerySymbolRuntime.Dynamic(CustomPropSymbol('membrane-orientation', 'is-transmembrane', Type.Bool),
    ctx => {
        const structure = ctx.currentStructure;
        const { x, y, z } = StructureProperties.atom;
        if (!structure.isAtomic) return 0;
        const membraneOrientation = MembraneOrientationProvider.get(structure).value;
        if (!membraneOrientation) return false;
        Vec3.set(pos, x(ctx.element), y(ctx.element), z(ctx.element));
        const { normal, p1, p2 } = membraneOrientation!;
        return isInMembranePlane(pos, normal, p1, p2);
    });

const isTransmembrane = StructureSelectionQuery('Residues embedded in membrane', MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'atom-test': isTransmembraneSymbol.symbol(),
            })
        ])
    ])
]), {
    description: 'Select residues that are embedded between the membrane layers.',
    category: StructureSelectionCategory.Residue,
    ensureCustomProperties: (ctx, structure) => {
        return MembraneOrientationProvider.attach(ctx, structure);
    }
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

        const colorTheme = AccessibleSurfaceAreaColorThemeProvider.name as any;
        return await PresetStructureRepresentations.auto.apply(ref, { ...params, globalThemeName: colorTheme }, plugin);
    }
});
