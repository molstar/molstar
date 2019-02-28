/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ResidueIndex, Model } from 'mol-model/structure';
import { BuiltInStructureRepresentationsName } from 'mol-repr/structure/registry';
import { BuiltInColorThemeName } from 'mol-theme/color';
import { AminoAcidNames } from 'mol-model/structure/model/types';
import { PluginContext } from 'mol-plugin/context';

export interface ModelInfo {
    hetResidues: { name: string, indices: ResidueIndex[] }[],
    assemblies: { id: string, details: string, isPreferred: boolean }[],
    preferredAssemblyId: string | undefined
}

export namespace ModelInfo {
    async function getPreferredAssembly(ctx: PluginContext, model: Model) {
        if (model.label.length <= 3) return void 0;
        try {
            const id = model.label.toLowerCase();
            const src = await ctx.runTask(ctx.fetch({ url: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/${id}` })) as string;
            const json = JSON.parse(src);
            const data = json && json[id];

            const assemblies = data[0] && data[0].assemblies;
            if (!assemblies || !assemblies.length) return void 0;

            for (const asm of assemblies) {
                if (asm.preferred) {
                    return asm.assembly_id;
                }
            }
            return void 0;
        } catch (e) {
            console.warn('getPreferredAssembly', e);
        }
    }

    export async function get(ctx: PluginContext, model: Model, checkPreferred: boolean): Promise<ModelInfo> {
        const { _rowCount: residueCount } = model.atomicHierarchy.residues;
        const { offsets: residueOffsets } = model.atomicHierarchy.residueAtomSegments;
        const chainIndex = model.atomicHierarchy.chainAtomSegments.index;
        // const resn = SP.residue.label_comp_id, entType = SP.entity.type;

        const pref = checkPreferred
            ? getPreferredAssembly(ctx, model)
            : void 0;

        const hetResidues: ModelInfo['hetResidues'] = [];
        const hetMap = new Map<string, ModelInfo['hetResidues'][0]>();

        for (let rI = 0 as ResidueIndex; rI < residueCount; rI++) {
            const comp_id = model.atomicHierarchy.residues.label_comp_id.value(rI);
            if (AminoAcidNames.has(comp_id)) continue;
            const mod_parent = model.properties.modifiedResidues.parentId.get(comp_id);
            if (mod_parent && AminoAcidNames.has(mod_parent)) continue;

            const cI = chainIndex[residueOffsets[rI]];
            const eI = model.atomicHierarchy.index.getEntityFromChain(cI);
            if (model.entities.data.type.value(eI) === 'water') continue;

            let lig = hetMap.get(comp_id);
            if (!lig) {
                lig = { name: comp_id, indices: [] };
                hetResidues.push(lig);
                hetMap.set(comp_id, lig);
            }
            lig.indices.push(rI);
        }

        const preferredAssemblyId = await pref;

        return {
            hetResidues: hetResidues,
            assemblies: model.symmetry.assemblies.map(a => ({ id: a.id, details: a.details, isPreferred: a.id === preferredAssemblyId })),
            preferredAssemblyId
        };
    }
}

export type SupportedFormats = 'cif' | 'pdb'
export interface LoadParams {
    url: string,
    format?: SupportedFormats,
    assemblyId?: string,
    representationStyle?: RepresentationStyle
}

export interface RepresentationStyle {
    sequence?: RepresentationStyle.Entry,
    hetGroups?: RepresentationStyle.Entry,
    water?: RepresentationStyle.Entry
}

export namespace RepresentationStyle {
    export type Entry = { kind?: BuiltInStructureRepresentationsName, coloring?: BuiltInColorThemeName }
}