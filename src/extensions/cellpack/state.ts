/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
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
import { Structure, StructureSymmetry, Unit, Model } from '../../mol-model/structure';
import { ModelSymmetry } from '../../mol-model-formats/structure/property/symmetry';
import { Vec3 } from '../../mol-math/linear-algebra';
import { StateTransformer } from '../../mol-state';
import { MBRepresentation, MBParams } from './representation';

export const DefaultCellPackBaseUrl = 'https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0';
export class CellPack extends PSO.Create<_CellPack>({ name: 'CellPack', typeClass: 'Object' }) { }

export { ParseCellPack };
type ParseCellPack = typeof ParseCellPack
const ParseCellPack = PluginStateTransform.BuiltIn({
    name: 'parse-cellpack',
    display: { name: 'Parse CellPack', description: 'Parse CellPack from JSON data' },
    from: PSO.Format.Json,
    to: CellPack,
    params: a => {
        return {
        };
    }
})({
    apply({a, params}) {
        return Task.create('Parse CellPack', async ctx => {
            const cell = a.data as Cell;
            let counter_id = 0;
            let fiber_counter_id = 0;
            let comp_counter = 0;
            const packings: CellPacking[] = [];
            const {compartments, cytoplasme } = cell;
            let iName = '';
            if(!cell.mapping_ids)cell.mapping_ids = {};
            if (cytoplasme) {
                packings.push({ name: 'Cytoplasme', location: 'cytoplasme', ingredients: cytoplasme.ingredients });
                for (iName in cytoplasme.ingredients){
                    if (cytoplasme.ingredients[iName].ingtype === 'fiber') {
                        cell.mapping_ids[-(fiber_counter_id + 1)] = [comp_counter, iName];
                        if (!cytoplasme.ingredients[iName].nbCurve) cytoplasme.ingredients[iName].nbCurve = 0;
                        fiber_counter_id++;
                    } else {
                        cell.mapping_ids[counter_id] = [comp_counter, iName];
                        if (!cytoplasme.ingredients[iName].results) {cytoplasme.ingredients[iName].results = [];}
                        counter_id++;
                    }
                }
                comp_counter++;
            }
            if (compartments) {
                for (const name in compartments) {
                    const { surface, interior } = compartments[name];
                    if (surface) {
                        packings.push({ name, location: 'surface', ingredients: surface.ingredients, geom: compartments[name].geom, geom_type: compartments[name].geom_type, mb: compartments[name].mb });
                        for (iName in surface.ingredients){
                            if (surface.ingredients[iName].ingtype === 'fiber') {
                                cell.mapping_ids[-(fiber_counter_id + 1)] = [comp_counter, iName];
                                if (!surface.ingredients[iName].nbCurve) surface.ingredients[iName].nbCurve = 0;
                                fiber_counter_id++;
                            } else {
                                cell.mapping_ids[counter_id] = [comp_counter, iName];
                                if (!surface.ingredients[iName].results) {surface.ingredients[iName].results = [];}
                                counter_id++;
                            }
                        }
                        comp_counter++;
                    }
                    if (interior) {
                        packings.push({ name, location: 'interior', ingredients: interior.ingredients });
                        for (iName in interior.ingredients){
                            if (interior.ingredients[iName].ingtype === 'fiber') {
                                cell.mapping_ids[-(fiber_counter_id + 1)] = [comp_counter, iName];
                                if (!interior.ingredients[iName].nbCurve) interior.ingredients[iName].nbCurve = 0;
                                fiber_counter_id++;
                            } else {
                                cell.mapping_ids[counter_id] = [comp_counter, iName];
                                if (!interior.ingredients[iName].results) {interior.ingredients[iName].results = [];}
                                counter_id++;
                            }
                        }
                        comp_counter++;
                    }
                }
            }
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
            const { structure, assets, colors } = await createStructureFromCellPack(plugin, packing, params.baseUrl, ingredientFiles).runInContext(ctx);
            await CellPackInfoProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset }, structure, {
                info: { packingsCount: a.data.packings.length, packingIndex: params.packing, colors }
            });
            (cache as any).assets = assets;
            return new PSO.Molecule.Structure(structure, { label: packing.name + '.' + packing.location });
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

export { StructureFromAssemblies };
type StructureFromAssemblies = typeof StructureFromAssemblies
const StructureFromAssemblies = PluginStateTransform.BuiltIn({
    name: 'Structure from all assemblies',
    display: { name: 'Structure from all assemblies' },
    from: PSO.Molecule.Model,
    to: PSO.Molecule.Structure,
    params: {
    }
})({
    canAutoUpdate({ newParams }) {
        return true;
    },
    apply({ a, params }) {
        return Task.create('Build Structure', async ctx => {
            // TODO: optimze
            // TODO: think of ways how to fast-track changes to this for animations
            const model = a.data;
            let initial_structure = Structure.ofModel(model);
            const structures: Structure[] = [];
            let structure: Structure = initial_structure;
            // the list of asambly *?
            const symmetry = ModelSymmetry.Provider.get(model);
            if (symmetry && symmetry.assemblies.length !== 0){
                for (const a of symmetry.assemblies) {
                    const s = await StructureSymmetry.buildAssembly(initial_structure, a.id).runInContext(ctx);
                    structures.push(s);
                }
                const builder = Structure.Builder({ label: 'Membrane' });
                let offsetInvariantId = 0;
                for (const s of structures) {
                    let maxInvariantId = 0;
                    for (const u of s.units) {
                        const invariantId = u.invariantId + offsetInvariantId;
                        if (u.invariantId > maxInvariantId) maxInvariantId = u.invariantId;
                        builder.addUnit(u.kind, u.model, u.conformation.operator, u.elements, Unit.Trait.None, invariantId);
                    }
                    offsetInvariantId += maxInvariantId + 1;
                }
                structure = builder.getStructure();
                for( let i = 0, il = structure.models.length; i < il; ++i) {
                    Model.TrajectoryInfo.set(structure.models[i], { size: il, index: i });
                }
            }
            return new PSO.Molecule.Structure(structure, { label: a.label, description: `${a.description}` });
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});


const CreateTransformer = StateTransformer.builderFactory('cellPACK');
export const CreateSphere = CreateTransformer({
    name: 'create-sphere',
    display: 'Sphere',
    from: PSO.Root, // or whatever data source
    to: PSO.Shape.Representation3D,
    params: {
        center: PD.Vec3(Vec3()),
        radius: PD.Numeric(1)
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Custom Sphere', async ctx => {
            const data = params;
            const repr = MBRepresentation({ webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },  () => (MBParams));
            await repr.createOrUpdate({ ...params, quality: 'custom', doubleSided:true }, data).runInContext(ctx);
            return new PSO.Shape.Representation3D({ repr, sourceData: a }, { label: `My Sphere` });
        });
    }
});