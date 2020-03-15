/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VisualQuality, VisualQualityOptions } from '../../../mol-geo/geometry/base';
import { Structure } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectRef } from '../../../mol-state';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StaticStructureComponentType } from '../../helpers/structure-component';
import { StructureSelectionQueries as Q } from '../../helpers/structure-selection-query';
import { PluginStateObject } from '../../objects';
import { RepresentationProviderTags, StructureRepresentationProvider } from './provider';
import { ColorTheme } from '../../../mol-theme/color';

export const CommonStructureRepresentationParams = {
    ignoreHydrogens: PD.Optional(PD.Boolean(false)),
    quality: PD.Optional(PD.Select<VisualQuality>('auto', VisualQualityOptions)),
    globalThemeName: PD.Optional(PD.Text<ColorTheme.BuiltIn>(''))
}
export type CommonStructureRepresentationParams = PD.ValuesFor<typeof CommonStructureRepresentationParams>

const auto = StructureRepresentationProvider({
    id: 'preset-structure-representation-auto',
    display: { name: 'Automatic', group: 'Preset' },
    params: () => CommonStructureRepresentationParams,
    apply(ctx, state, structureCell, params, plugin) {
        const structure = structureCell.obj!.data;
        const size = Structure.getSize(structure)

        switch (size) {
            case Structure.Size.Gigantic:
            case Structure.Size.Huge:
                return coarseSurface.apply(ctx, state, structureCell, params, plugin);
            case Structure.Size.Large:
                return polymerCartoon.apply(ctx, state, structureCell, params, plugin);
            case Structure.Size.Medium:
                return polymerAndLigand.apply(ctx, state, structureCell, params, plugin);
            case Structure.Size.Small:
                return atomicDetail.apply(ctx, state, structureCell, params, plugin);
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

const polymerAndLigand = StructureRepresentationProvider({
    id: 'preset-structure-representation-polymer-and-ligand',
    display: { name: 'Polymer & Ligand', group: 'Preset' },
    params: () => CommonStructureRepresentationParams,
    async apply(ctx, state, structureCell, params, plugin) {
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

        await state.updateTree(update, { revertOnError: false }).runInContext(ctx);
        return { components, representations };
    }
});

const proteinAndNucleic = StructureRepresentationProvider({
    id: 'preset-structure-representation-protein-and-nucleic',
    display: { name: 'Protein & Nucleic', group: 'Preset' },
    params: () => CommonStructureRepresentationParams,
    async apply(ctx, state, structureCell, params, plugin) {
        const components = {
            protein: await presetSelectionComponent(plugin, structureCell, 'protein'),
            nucleic: await presetSelectionComponent(plugin, structureCell, 'nucleic'),
        };

        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            protein: builder.buildRepresentation(update, components.protein, { type: 'cartoon', typeParams, color }),
            nucleic: builder.buildRepresentation(update, components.nucleic, { type: 'gaussian-surface', typeParams, color })
        };

        await state.updateTree(update, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const coarseSurface = StructureRepresentationProvider({
    id: 'preset-structure-representation-coarse-surface',
    display: { name: 'Coarse Surface', group: 'Preset' },
    params: () => CommonStructureRepresentationParams,
    async apply(ctx, state, structureCell, params, plugin) {
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

        await state.updateTree(update, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const polymerCartoon = StructureRepresentationProvider({
    id: 'preset-structure-representation-polymer-cartoon',
    display: { name: 'Polymer Cartoon', group: 'Preset' },
    params: () => CommonStructureRepresentationParams,
    async apply(ctx, state, structureCell, params, plugin) {
        const components = {
            polymer: await presetSelectionComponent(plugin, structureCell, 'polymer'),
        };

        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            polymer: builder.buildRepresentation(update, components.polymer, { type: 'cartoon', typeParams, color })
        };

        await state.updateTree(update, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const atomicDetail = StructureRepresentationProvider({
    id: 'preset-structure-representation-atomic-detail',
    display: { name: 'Atomic Detail', group: 'Preset' },
    params: () => CommonStructureRepresentationParams,
    async apply(ctx, state, structureCell, params, plugin) {

        const components = {
            all: await presetSelectionComponent(plugin, structureCell, 'all'),
        };

        const { update, builder, typeParams, color } = reprBuilder(plugin, params);
        const representations = {
            all: builder.buildRepresentation(update, components.all, { type: 'ball-and-stick', typeParams, color })
        };

        await state.updateTree(update, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

export function presetStaticComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, type: StaticStructureComponentType) {
    return plugin.builders.structure.tryCreateStaticComponent({
        structure,
        type,
        key: `static-${type}`,
        tags: [RepresentationProviderTags.Component]
    });
}

export function presetSelectionComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, query: keyof typeof Q) {
    return plugin.builders.structure.tryCreateQueryComponent({ 
        structure,
        query: Q[query],
        key: `selection-${query}`,
        tags: [RepresentationProviderTags.Component]
    });
}

export const PresetStructureReprentations = {
    auto,
    atomicDetail,
    polymerCartoon,
    polymerAndLigand,
    proteinAndNucleic,
    coarseSurface
};
export type PresetStructureReprentations = typeof PresetStructureReprentations;