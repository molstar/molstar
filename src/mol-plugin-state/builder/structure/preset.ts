/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransforms } from '../../transforms';
import { StructureRepresentation3DHelpers } from '../../transforms/representation';
import { StructureSelectionQueries as Q } from '../../../mol-plugin/util/structure-selection-helper';
import { BuiltInStructureRepresentations } from '../../../mol-repr/structure/registry';
import { StructureRepresentationProvider, RepresentationProviderTags } from './provider';
import { StateObjectRef } from '../../../mol-state';
import { PluginStateObject } from '../../objects';
import { StaticStructureComponentType } from '../../helpers/structure-component';
import { PluginContext } from '../../../mol-plugin/context';
import { Structure } from '../../../mol-model/structure';

const auto = StructureRepresentationProvider({
    id: 'preset-structure-representation-auto',
    display: { name: 'Automatic', group: 'Preset' },
    apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj!.data;
        const size = Structure.getSize(structure)

        switch (size) {
            case Structure.Size.Gigantic:
            case Structure.Size.Huge:
                return coarseSurface.apply(ctx, state, structureCell, void 0, plugin);
            case Structure.Size.Large:
                return polymerCartoon.apply(ctx, state, structureCell, void 0, plugin);
            case Structure.Size.Medium:
                return polymerAndLigand.apply(ctx, state, structureCell, void 0, plugin);
            case Structure.Size.Small:
                return atomicDetail.apply(ctx, state, structureCell, void 0, plugin);
        }
    }
});

const polymerAndLigand = StructureRepresentationProvider({
    id: 'preset-structure-representation-polymer-and-ligand',
    display: { name: 'Polymer & Ligand', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj?.data!;
        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const components = {
            polymer: await staticComponent(plugin, structureCell, 'polymer'),
            ligand: await staticComponent(plugin, structureCell, 'ligand'),
            nonStandard: await staticComponent(plugin, structureCell, 'non-standard'),
            branched: await staticComponent(plugin, structureCell, 'branched'),
            water: await staticComponent(plugin, structureCell, 'water'),
            coarse: await staticComponent(plugin, structureCell, 'coarse')
        };

        const builder = state.build();
        const representations = {
            polymer: components.polymer && builder
                .to(components.polymer)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure)).selector,
            ligand: components.ligand && builder
                .to(components.ligand)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure)).selector,
            nonStandard: components.nonStandard && builder
                .to(components.nonStandard)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'ball-and-stick', 'polymer-id', structure, void 0)).selector,
            branched: components.branched && {
                ballAndStick: builder.to(components.branched).applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.15 })).selector,
                snfg3d: builder.to(components.branched).applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'carbohydrate', structure)).selector
            },
            water: components.water && builder
                .to(components.water)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.51 })).selector,
            coarse: components.coarse && builder
                .to(components.coarse)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'spacefill', 'polymer-id', structure, {}))
        };

        await state.updateTree(builder, { revertOnError: false }).runInContext(ctx);
        return { components, representations };
    }
});

const proteinAndNucleic = StructureRepresentationProvider({
    id: 'preset-structure-representation-protein-and-nucleic',
    display: { name: 'Protein & Nucleic', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj!.data;
        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const components = {
            protein: await selectionComponent(plugin, structureCell, 'protein'),
            nucleic: await selectionComponent(plugin, structureCell, 'nucleic'),
        };

        const builder = state.build();
        const representations = {
            protein: components.protein && builder
                .to(components.protein)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure)).selector,
            nucleic: components.nucleic && builder
                .to(components.nucleic)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'gaussian-surface', structure)).selector,
        };

        await state.updateTree(builder, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const coarseSurface = StructureRepresentationProvider({
    id: 'preset-structure-representation-coarse-surface',
    display: { name: 'Coarse Surface', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
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
            components.trace = await selectionComponent(plugin, structureCell, 'trace')
        } else if(size === Structure.Size.Huge) {
            Object.assign(gaussianProps, {
                smoothness: 0.5,
            })
            components.trace = await selectionComponent(plugin, structureCell, 'polymer')
        } else {
            components.trace = await selectionComponent(plugin, structureCell, 'polymer')
        }

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: [
                BuiltInStructureRepresentations['gaussian-surface'],
                () => gaussianProps
            ]
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const builder = state.build();
        const representations = {
            trace: components.trace && builder
                .to(components.trace)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params).selector,
        };

        await state.updateTree(builder, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const polymerCartoon = StructureRepresentationProvider({
    id: 'preset-structure-representation-polymer-cartoon',
    display: { name: 'Polymer Cartoon', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj!.data;
        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const components = {
            polymer: await selectionComponent(plugin, structureCell, 'polymer'),
        };

        const builder = state.build();
        const representations = {
            polymer: components.polymer && builder
                .to(components.polymer)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure)).selector,
        };

        await state.updateTree(builder, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const atomicDetail = StructureRepresentationProvider({
    id: 'preset-structure-representation-atomic-detail',
    display: { name: 'Atomic Detail', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj!.data;
        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const components = {
            all: await selectionComponent(plugin, structureCell, 'all'),
        };

        const builder = state.build();
        const representations = {
            all: components.all && builder
                .to(components.all)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure)).selector,
        };

        await state.updateTree(builder, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

function staticComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, type: StaticStructureComponentType) {
    return plugin.builders.structure.tryCreateComponent(structure, {
        type: { name: 'static', params: type },
        nullIfEmpty: true,
        label: ''
    }, `static-${type}`);
}

function selectionComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, query: keyof typeof Q) {
    return plugin.builders.structure.tryCreateComponent(structure, {
        type: { name: 'expression', params: Q[query].expression },
        nullIfEmpty: true,
        label: Q[query].label
    }, `selection-${query}`);
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