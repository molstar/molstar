/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject as PSO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Task } from '../../mol-task';
import { CellPack as _CellPack, Cell, CellPacking } from './data';
import { createStructureFromCellPack } from './model';
import { IngredientFiles } from './util';
import { Asset } from '../../mol-util/assets';
import { PluginContext } from '../../mol-plugin/context';
import { CellPackInfoProvider } from './property';

export const DefaultCellPackBaseUrl = 'https://mesoscope.scripps.edu/data/cellPACK_data/cellPACK_database_1.1.0/';

export class CellPack extends PSO.Create<_CellPack>({ name: 'CellPack', typeClass: 'Object' }) { }

export { ParseCellPack };
type ParseCellPack = typeof ParseCellPack
const ParseCellPack = PluginStateTransform.BuiltIn({
    name: 'parse-cellpack',
    display: { name: 'Parse CellPack', description: 'Parse CellPack from JSON data' },
    from: PSO.Format.Json,
    to: CellPack
})({
    apply({ a }) {
        return Task.create('Parse CellPack', async ctx => {
            const cell = a.data as Cell;

            const packings: CellPacking[] = [];
            const { compartments, cytoplasme } = cell;
            if (compartments) {
                for (const name in compartments) {
                    const { surface, interior } = compartments[name];
                    if (surface) packings.push({ name, location: 'surface', ingredients: surface.ingredients });
                    if (interior) packings.push({ name, location: 'interior', ingredients: interior.ingredients });
                }
            }
            if (cytoplasme) packings.push({ name: 'Cytoplasme', location: 'cytoplasme', ingredients: cytoplasme.ingredients });

            return new CellPack({ cell, packings });
        });
    }
});

export { StructureFromCellpack };
type StructureFromCellpack = typeof ParseCellPack
const StructureFromCellpack = PluginStateTransform.BuiltIn({
    name: 'structure-from-cellpack',
    display: { name: 'Structure from CellPack', description: 'Create Structure from CellPack Packing' },
    from: CellPack,
    to: PSO.Molecule.Structure,
    params: a => {
        const options = a ? a.data.packings.map((d, i) => [i, d.name] as const) : [];
        return {
            packing: PD.Select(0, options),
            baseUrl: PD.Text(DefaultCellPackBaseUrl),
            ingredientFiles: PD.FileList({ accept: '.cif,.bcif,.pdb' })
        };
    }
})({
    apply({ a, params, cache }, plugin: PluginContext) {
        return Task.create('Structure from CellPack', async ctx => {
            const packing = a.data.packings[params.packing];
            const ingredientFiles: IngredientFiles = {};
            if (params.ingredientFiles !== null) {
                for (const file of params.ingredientFiles) {
                    ingredientFiles[file.name] = file;
                }
            }
            const { structure, assets } = await createStructureFromCellPack(plugin.managers.asset, packing, params.baseUrl, ingredientFiles).runInContext(ctx);

            await CellPackInfoProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset }, structure, {
                info: { packingsCount: a.data.packings.length, packingIndex: params.packing }
            });

            (cache as any).assets = assets;
            return new PSO.Molecule.Structure(structure, { label: packing.name });
        });
    },
    dispose({ b, cache }) {
        const assets = (cache as any).assets as Asset.Wrapper[];
        if(assets) {
            for (const a of assets) a.dispose();
        }

        if (b) {
            b.data.customPropertyDescriptors.dispose();
            for (const m of b.data.models) {
                m.customProperties.dispose();
            }
        }
    }
});