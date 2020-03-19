/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PresetProvider } from '../preset-provider';
import { PluginStateObject } from '../../objects';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualQuality, VisualQualityOptions } from '../../../mol-geo/geometry/base';
import { ColorTheme } from '../../../mol-theme/color';
import { Structure } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectRef } from '../../../mol-state';
import { StaticStructureComponentType } from '../../helpers/structure-component';
import { StructureSelectionQueries as Q } from '../../helpers/structure-selection-query';

export interface StructureRepresentationPresetProvider<P = any, S = {}> extends PresetProvider<PluginStateObject.Molecule.Structure, P, S> { }
export namespace StructureRepresentationPresetProvider {
    export type Params<P extends StructureRepresentationPresetProvider> = P extends StructureRepresentationPresetProvider<infer T> ? T : never;
    export type State<P extends StructureRepresentationPresetProvider> = P extends StructureRepresentationPresetProvider<infer _, infer S> ? S : never;
}
export function StructureRepresentationPresetProvider<P, S>(repr: StructureRepresentationPresetProvider<P, S>) { return repr; }

export const CommonStructureRepresentationPresetParams = {
    ignoreHydrogens: PD.Optional(PD.Boolean(false)),
    quality: PD.Optional(PD.Select<VisualQuality>('auto', VisualQualityOptions)),
    globalThemeName: PD.Optional(PD.Text<ColorTheme.BuiltIn>(''))
}
export type CommonStructureRepresentationParams = PD.ValuesFor<typeof CommonStructureRepresentationPresetParams>

const auto = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-auto',
    display: { name: 'Automatic', group: 'Preset' },
    params: () => CommonStructureRepresentationPresetParams,
    apply(ref, params, plugin) {
        const structure = StateObjectRef.resolveAndCheck(plugin.state.data, ref)?.obj?.data;
        if (!structure) return { };
        const size = Structure.getSize(structure)

        switch (size) {
            case Structure.Size.Gigantic:
            case Structure.Size.Huge:
                return coarseSurface.apply(ref, params, plugin);
            case Structure.Size.Large:
                return polymerCartoon.apply(ref, params, plugin);
            case Structure.Size.Medium:
                return polymerAndLigand.apply(ref, params, plugin);
            case Structure.Size.Small:
                return atomicDetail.apply(ref, params, plugin);
        }
    }
});

function reprBuilder(plugin: PluginContext, params: CommonStructureRepresentationParams) {
    const update = plugin.state.data.build();
    const builder = plugin.builders.structure.representation;
    const typeParams = {
        quality: plugin.managers.structure.component.state.options.visualQuality,
        ignoreHydrogens: !plugin.managers.structure.component.state.options.showHydrogens,
    };
    if (params.quality && params.quality !== 'auto') typeParams.quality = params.quality;
    if (params.ignoreHydrogens !== void 0) typeParams.ignoreHydrogens = !!params.ignoreHydrogens;
    const color: ColorTheme.BuiltIn | undefined = params.globalThemeName ? params.globalThemeName : void 0;

    return { update, builder, color, typeParams };
}

const empty = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-empty',
    display: { name: 'Empty', group: 'Preset' },
    async apply(ref, params, plugin) {
        return { };
    }
});

const polymerAndLigand = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-polymer-and-ligand',
    display: { name: 'Polymer & Ligand', group: 'Preset' },
    params: () => CommonStructureRepresentationPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const components = {
            polymer: await presetStaticComponent(plugin, structureCell, 'polymer'),
            ligand: await presetStaticComponent(plugin, structureCell, 'ligand'),
            nonStandard: await presetStaticComponent(plugin, structureCell, 'non-standard'),
            branched: await presetStaticComponent(plugin, structureCell, 'branched'),
            water: await presetStaticComponent(plugin, structureCell, 'water'),
            coarse: await presetStaticComponent(plugin, structureCell, 'coarse')
        };

        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            polymer: builder.buildRepresentation(update, components.polymer, { type: 'cartoon', typeParams, color }),
            ligand: builder.buildRepresentation(update, components.ligand, { type: 'ball-and-stick', typeParams, color }),
            nonStandard: builder.buildRepresentation(update, components.nonStandard, { type: 'ball-and-stick', typeParams, color: color || 'polymer-id' }),
            branched: components.branched && {
                ballAndStick: builder.buildRepresentation(update, components.branched, { type: 'ball-and-stick', typeParams: { ...typeParams, alpha: 0.15 }, color }),
                snfg3d: builder.buildRepresentation(update, components.branched, { type: 'carbohydrate', typeParams, color }),
            },
            water: builder.buildRepresentation(update, components.water, { type: 'ball-and-stick', typeParams: { ...typeParams, alpha: 0.51 }, color }),
            coarse: builder.buildRepresentation(update, components.coarse, { type: 'spacefill', typeParams, color: color || 'polymer-id' })
        };

        await plugin.updateDataState(update, { revertOnError: false });
        return { components, representations };
    }
});

const proteinAndNucleic = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-protein-and-nucleic',
    display: { name: 'Protein & Nucleic', group: 'Preset' },
    params: () => CommonStructureRepresentationPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const components = {
            protein: await presetSelectionComponent(plugin, structureCell, 'protein'),
            nucleic: await presetSelectionComponent(plugin, structureCell, 'nucleic'),
        };

        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            protein: builder.buildRepresentation(update, components.protein, { type: 'cartoon', typeParams, color }),
            nucleic: builder.buildRepresentation(update, components.nucleic, { type: 'gaussian-surface', typeParams, color })
        };

        await plugin.updateDataState(update, { revertOnError: true });
        return { components, representations };
    }
});

const coarseSurface = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-coarse-surface',
    display: { name: 'Coarse Surface', group: 'Preset' },
    params: () => CommonStructureRepresentationPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const structure = structureCell.obj!.data;
        const size = Structure.getSize(structure)

        const gaussianProps = Object.create(null);
        const components = Object.create(null);

        if (size === Structure.Size.Gigantic) {
            Object.assign(gaussianProps, {
                radiusOffset: 1,
                smoothness: 0.5,
                visuals: ['structure-gaussian-surface-mesh']
            })
            components.trace = await presetSelectionComponent(plugin, structureCell, 'trace')
        } else if(size === Structure.Size.Huge) {
            Object.assign(gaussianProps, {
                smoothness: 0.5,
            })
            components.trace = await presetSelectionComponent(plugin, structureCell, 'polymer')
        } else {
            components.trace = await presetSelectionComponent(plugin, structureCell, 'polymer')
        }


        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            trace: builder.buildRepresentation(update, components.trace, { type: 'gaussian-surface', typeParams: { ...typeParams, ...gaussianProps }, color })
        };

        await plugin.updateDataState(update, { revertOnError: true });
        return { components, representations };
    }
});

const polymerCartoon = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-polymer-cartoon',
    display: { name: 'Polymer Cartoon', group: 'Preset' },
    params: () => CommonStructureRepresentationPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const components = {
            polymer: await presetSelectionComponent(plugin, structureCell, 'polymer'),
        };

        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            polymer: builder.buildRepresentation(update, components.polymer, { type: 'cartoon', typeParams, color })
        };

        await plugin.updateDataState(update, { revertOnError: true });
        return { components, representations };
    }
});

const atomicDetail = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-atomic-detail',
    display: { name: 'Atomic Detail', group: 'Preset' },
    params: () => CommonStructureRepresentationPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const components = {
            all: await presetSelectionComponent(plugin, structureCell, 'all'),
        };

        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            all: builder.buildRepresentation(update, components.all, { type: 'ball-and-stick', typeParams, color })
        };

        await plugin.updateDataState(update, { revertOnError: true });
        return { components, representations };
    }
});

export function presetStaticComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, type: StaticStructureComponentType) {
    return plugin.builders.structure.tryCreateStaticComponent({
        structure,
        type,
        key: `static-${type}`
    });
}

export function presetSelectionComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, query: keyof typeof Q) {
    return plugin.builders.structure.tryCreateQueryComponent({
        structure,
        query: Q[query],
        key: `selection-${query}`
    });
}

export const PresetStructureReprentations = {
    empty,
    auto,
    'atomic-detail': atomicDetail,
    'polymer-cartoon': polymerCartoon,
    'polymer-and-ligand': polymerAndLigand,
    'protein-and-nucleic': proteinAndNucleic,
    'coarse-surface': coarseSurface
};
export type PresetStructureReprentations = typeof PresetStructureReprentations;