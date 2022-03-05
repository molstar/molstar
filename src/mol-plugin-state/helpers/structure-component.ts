/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Expression } from '../../mol-script/language/expression';
import { MolScriptBuilder } from '../../mol-script/language/builder';
import { StructureElement, Structure, StructureSelection as Sel, StructureQuery, Queries, QueryContext, Model } from '../../mol-model/structure';
import { StructureQueryHelper } from './structure-query';
import { PluginStateObject as SO } from '../objects';
import { StructureSelectionQueries } from './structure-selection-query';
import { StateTransformer, StateObject } from '../../mol-state';
import { Script } from '../../mol-script/script';
import { assertUnreachable } from '../../mol-util/type-helpers';

export const StaticStructureComponentTypes = [
    'all',

    'polymer',

    'protein',
    'nucleic',
    'water',
    'ion',
    'lipid',

    'branched',
    'ligand',
    'non-standard',

    'coarse'
] as const;

export type StaticStructureComponentType = (typeof StaticStructureComponentTypes)[number]

export const StructureComponentParams = () => ({
    type: PD.MappedStatic('static', {
        static: PD.Text<StaticStructureComponentType>('polymer'),
        expression: PD.Value<Expression>(MolScriptBuilder.struct.generator.all),
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty),
        script: PD.Script({ language: 'mol-script', expression: '(sel.atom.all)' }),
    }, { isHidden: true }),
    nullIfEmpty: PD.Optional(PD.Boolean(true, { isHidden: true })),
    label: PD.Text('', { isHidden: true })
});
export type StructureComponentParams = PD.ValuesFor<ReturnType<typeof StructureComponentParams>>

export function createStructureComponent(a: Structure, params: StructureComponentParams, cache: { source: Structure, entry?: StructureQueryHelper.CacheEntry }) {
    cache.source = a;

    let component: Structure = Structure.Empty;
    let label: string | undefined = void 0;
    switch (params.type.name) {
        case 'static': {
            let query: StructureQuery;
            switch (params.type.params) {
                case 'all': query = StructureSelectionQueries.all.query; label = 'All'; break;

                case 'polymer': query = StructureSelectionQueries.polymer.query; label = 'Polymer'; break;

                case 'protein': query = StructureSelectionQueries.protein.query; label = 'Protein'; break;
                case 'nucleic': query = StructureSelectionQueries.nucleic.query; label = 'Nucleic'; break;
                case 'water': query = Queries.internal.water(); label = 'Water'; break;
                case 'ion': query = StructureSelectionQueries.ion.query; label = 'Ion'; break;
                case 'lipid': query = StructureSelectionQueries.lipid.query; label = 'Lipid'; break;

                case 'branched': query = StructureSelectionQueries.branchedPlusConnected.query; label = 'Branched'; break;
                case 'ligand': query = StructureSelectionQueries.ligandPlusConnected.query; label = 'Ligand'; break;

                case 'non-standard': query = StructureSelectionQueries.nonStandardPolymer.query; label = 'Non-standard'; break;

                case 'coarse': query = StructureSelectionQueries.coarse.query; label = 'Coarse'; break;

                default: assertUnreachable(params.type);
            }
            const result = query(new QueryContext(a));
            component = Sel.unionStructure(result);
            break;
        }
        case 'script':
        case 'expression': {
            const { selection, entry } = StructureQueryHelper.createAndRun(a, params.type.params);
            cache.entry = entry;
            component = Sel.unionStructure(selection);
            break;
        }
        case 'bundle': {
            if (params.type.params.hash !== a.hashCode) break;
            component = StructureElement.Bundle.toStructure(params.type.params, a);
            break;
        }
    }

    if (params.nullIfEmpty && component.elementCount === 0) return StateObject.Null;

    const props = { label: `${params.label || label || 'Component'}`, description: Structure.elementDescription(component) };
    return new SO.Molecule.Structure(component, props);
}

export function updateStructureComponent(a: Structure, b: SO.Molecule.Structure, oldParams: StructureComponentParams, newParams: StructureComponentParams, cache: { source: Structure, entry?: StructureQueryHelper.CacheEntry }) {
    if (oldParams.type.name !== newParams.type.name) return StateTransformer.UpdateResult.Recreate;

    let updated = false;

    switch (newParams.type.name) {
        case 'static': {
            if (oldParams.type.params !== newParams.type.params) {
                return StateTransformer.UpdateResult.Recreate;
            }
            if (!Structure.areEquivalent(a, cache.source)) {
                return StateTransformer.UpdateResult.Recreate;
            }
            if (b.data.model === a.model) return StateTransformer.UpdateResult.Unchanged;
            if (!Model.areHierarchiesEqual(a.model, b.data.model)) return StateTransformer.UpdateResult.Recreate;

            b.data = b.data.remapModel(a.model);
            return StateTransformer.UpdateResult.Updated;
        }
        case 'script':
            if (!Script.areEqual(oldParams.type.params as Script, newParams.type.params)) {
                return StateTransformer.UpdateResult.Recreate;
            }
        case 'expression': {
            if ((oldParams.type.params as Expression) !== newParams.type.params) {
                return StateTransformer.UpdateResult.Recreate;
            }

            if (a === cache.source) break;

            const entry = (cache as { entry: StructureQueryHelper.CacheEntry }).entry;

            const selection = StructureQueryHelper.updateStructure(entry, a);
            cache.source = a;
            b.data = Sel.unionStructure(selection);
            StructureQueryHelper.updateStructureObject(b, selection, newParams.label);
            updated = true;
            break;
        }
        case 'bundle': {
            if (a === cache.source && StructureElement.Bundle.areEqual(oldParams.type.params as StructureElement.Bundle, newParams.type.params)) {
                break;
            }

            cache.source = a;
            if (newParams.type.params.hash !== a.hashCode) {
                updated = b.data.elementCount !== 0;
                b.data = b.data.elementCount === 0 ? b.data : Structure.Empty;
            } else {
                updated = true;
                b.data = StructureElement.Bundle.toStructure(newParams.type.params, a);
            }
            break;
        }
    }

    if (updated) {
        if (newParams.nullIfEmpty && b.data.elementCount === 0) return StateTransformer.UpdateResult.Null;

        b.description = Structure.elementDescription(b.data);
    }

    if (oldParams.label !== newParams.label) {
        updated = true;
        b.label = `${newParams.label || b.label}`;
    }

    return updated ? StateTransformer.UpdateResult.Updated : StateTransformer.UpdateResult.Unchanged;
}