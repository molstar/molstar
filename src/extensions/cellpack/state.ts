/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from '../../mol-data/int';
import { flipByteOrder, IsNativeEndianLittle } from '../../mol-io/common/binary';
import { SymmetryOperator } from '../../mol-math/geometry';
import { Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { PluginStateObject as PSO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CellPack as _CellPack, Cell, CellPacking } from './data';
import { createStructureFromCellPack } from './model';
import { getFloatValue, IngredientFiles } from './util';
import { Asset } from '../../mol-util/assets';
import { PluginContext } from '../../mol-plugin/context';
import { CellPackInfoProvider } from './property';
import { ModelSymmetry } from '../../mol-model-formats/structure/property/symmetry';
import { ElementIndex, EntityIndex, Model, Structure, Unit } from '../../mol-model/structure';
import { Assembly, Symmetry } from '../../mol-model/structure/model/properties/symmetry';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { MBParams, MBRepresentation } from './representation';

export const DefaultCellPackBaseUrl = 'https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0';
export class CellPack extends PSO.Create<_CellPack>({ name: 'CellPack', typeClass: 'Object' }) { }

export { ParseCellPack };
export { StructureFromCellpack };
export { CellpackAssembly };
export { EntityStructure };
type ParseCellPack = typeof ParseCellPack
const ParseCellPack = PluginStateTransform.BuiltIn({
    name: 'parse-cellpack',
    display: { name: 'Parse CellPack', description: 'Parse CellPack from JSON data' },
    from: PSO.Format.Json,
    to: CellPack,
    params: a => {
        return {
            resultsFile: PD.File({ accept: '.bin' }),
            baseUrl: PD.Text(DefaultCellPackBaseUrl)
        };
    }
})({
    apply({ a, params, cache }, plugin: PluginContext) {
        return Task.create('Parse CellPack', async ctx => {
            const cell = a.data as Cell;
            let counter_id = 0;
            let fiber_counter_id = 0;
            let comp_counter = 0;
            const packings: CellPacking[] = [];
            const { compartments, cytoplasme } = cell;
            if (!cell.mapping_ids) cell.mapping_ids = {};
            if (cytoplasme) {
                packings.push({ name: 'Cytoplasme', location: 'cytoplasme', ingredients: cytoplasme.ingredients });
                for (const iName in cytoplasme.ingredients) {
                    if (cytoplasme.ingredients[iName].ingtype === 'fiber') {
                        cell.mapping_ids[-(fiber_counter_id + 1)] = [comp_counter, iName];
                        if (!cytoplasme.ingredients[iName].nbCurve) cytoplasme.ingredients[iName].nbCurve = 0;
                        fiber_counter_id++;
                    } else {
                        cell.mapping_ids[counter_id] = [comp_counter, iName];
                        if (!cytoplasme.ingredients[iName].results) { cytoplasme.ingredients[iName].results = []; }
                        counter_id++;
                    }
                }
                comp_counter++;
            }
            if (compartments) {
                for (const name in compartments) {
                    const { surface, interior } = compartments[name];
                    let filename = '';
                    if (compartments[name].geom_type === 'file') {
                        filename = (compartments[name].geom) ? compartments[name].geom as string : '';
                    }
                    const compartment = { filename: filename, geom_type: compartments[name].geom_type, compartment_primitives: compartments[name].mb };
                    if (surface) {
                        packings.push({ name, location: 'surface', ingredients: surface.ingredients, compartment: compartment });
                        for (const iName in surface.ingredients) {
                            if (surface.ingredients[iName].ingtype === 'fiber') {
                                cell.mapping_ids[-(fiber_counter_id + 1)] = [comp_counter, iName];
                                if (!surface.ingredients[iName].nbCurve) surface.ingredients[iName].nbCurve = 0;
                                fiber_counter_id++;
                            } else {
                                cell.mapping_ids[counter_id] = [comp_counter, iName];
                                if (!surface.ingredients[iName].results) { surface.ingredients[iName].results = []; }
                                counter_id++;
                            }
                        }
                        comp_counter++;
                    }
                    if (interior) {
                        if (!surface) packings.push({ name, location: 'interior', ingredients: interior.ingredients, compartment: compartment });
                        else packings.push({ name, location: 'interior', ingredients: interior.ingredients });
                        for (const iName in interior.ingredients) {
                            if (interior.ingredients[iName].ingtype === 'fiber') {
                                cell.mapping_ids[-(fiber_counter_id + 1)] = [comp_counter, iName];
                                if (!interior.ingredients[iName].nbCurve) interior.ingredients[iName].nbCurve = 0;
                                fiber_counter_id++;
                            } else {
                                cell.mapping_ids[counter_id] = [comp_counter, iName];
                                if (!interior.ingredients[iName].results) { interior.ingredients[iName].results = []; }
                                counter_id++;
                            }
                        }
                        comp_counter++;
                    }
                }
            }
            const { options } = cell;
            let resultsAsset: Asset.Wrapper<'binary'> | undefined;
            if (params.resultsFile) {
                resultsAsset = await plugin.runTask(plugin.managers.asset.resolve(params.resultsFile, 'binary', true));
            } else if (options?.resultfile) {
                const url = `${params.baseUrl}/results/${options.resultfile}`;
                resultsAsset = await plugin.runTask(plugin.managers.asset.resolve(Asset.getUrlAsset(plugin.managers.asset, url), 'binary', true));
            }
            if (resultsAsset) {
                (cache as any).asset = resultsAsset;
                const results = resultsAsset.data;
                // flip the byte order if needed
                const buffer = IsNativeEndianLittle ? results.buffer : flipByteOrder(results, 4);
                const numbers = new DataView(buffer);
                const ninst = getFloatValue(numbers, 0);
                const npoints = getFloatValue(numbers, 4);
                const ncurve = getFloatValue(numbers, 8);

                let offset = 12;

                if (ninst !== 0) {
                    const pos = new Float32Array(buffer, offset, ninst * 4);
                    offset += ninst * 4 * 4;
                    const quat = new Float32Array(buffer, offset, ninst * 4);
                    offset += ninst * 4 * 4;

                    for (let i = 0; i < ninst; i++) {
                        const x: number = pos[i * 4 + 0];
                        const y: number = pos[i * 4 + 1];
                        const z: number = pos[i * 4 + 2];
                        const ingr_id = pos[i * 4 + 3] as number;
                        const pid = cell.mapping_ids![ingr_id];
                        if (!packings[pid[0]].ingredients[pid[1]].results) {
                            packings[pid[0]].ingredients[pid[1]].results = [];
                        }
                        packings[pid[0]].ingredients[pid[1]].results.push([Vec3.create(x, y, z),
                            Quat.create(quat[i * 4 + 0], quat[i * 4 + 1], quat[i * 4 + 2], quat[i * 4 + 3])]);
                    }
                }

                if (npoints !== 0) {
                    const ctr_pos = new Float32Array(buffer, offset, npoints * 4);
                    offset += npoints * 4 * 4;
                    offset += npoints * 4 * 4;
                    const ctr_info = new Float32Array(buffer, offset, npoints * 4);
                    offset += npoints * 4 * 4;
                    const curve_ids = new Float32Array(buffer, offset, ncurve * 4);
                    offset += ncurve * 4 * 4;

                    let counter = 0;
                    let ctr_points: Vec3[] = [];
                    let prev_ctype = 0;
                    let prev_cid = 0;

                    for (let i = 0; i < npoints; i++) {
                        const x: number = -ctr_pos[i * 4 + 0];
                        const y: number = ctr_pos[i * 4 + 1];
                        const z: number = ctr_pos[i * 4 + 2];
                        const cid: number = ctr_info[i * 4 + 0]; // curve id
                        const ctype: number = curve_ids[cid * 4 + 0]; // curve type
                        // cid  148 165 -1 0
                        // console.log("cid ",cid,ctype,prev_cid,prev_ctype);//165,148
                        if (prev_ctype !== ctype) {
                            const pid = cell.mapping_ids![-prev_ctype - 1];
                            const cname = `curve${counter}`;
                            packings[pid[0]].ingredients[pid[1]].nbCurve = counter + 1;
                            packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
                            ctr_points = [];
                            counter = 0;
                        } else if (prev_cid !== cid) {
                            ctr_points = [];
                            const pid = cell.mapping_ids![-prev_ctype - 1];
                            const cname = `curve${counter}`;
                            packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
                            counter += 1;
                        }
                        ctr_points.push(Vec3.create(x, y, z));
                        prev_ctype = ctype;
                        prev_cid = cid;
                    }

                    // do the last one
                    const pid = cell.mapping_ids![-prev_ctype - 1];
                    const cname = `curve${counter}`;
                    packings[pid[0]].ingredients[pid[1]].nbCurve = counter + 1;
                    packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
                }
            }
            return new CellPack({ cell, packings });
        });
    },
    dispose({ cache }) {
        ((cache as any)?.asset as Asset.Wrapper | undefined)?.dispose();
    },
});

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
        if (assets) {
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

function createModelChainMap(model: Model) {
    const builder = new Structure.StructureBuilder();
    const units = new Map<string, Unit>();

    const { label_asym_id, _rowCount } = model.atomicHierarchy.chains;
    const { offsets } = model.atomicHierarchy.chainAtomSegments;

    for (let i = 0; i < _rowCount; i++) {
        const elements = SortedArray.ofBounds(offsets[i] as ElementIndex, offsets[i + 1] as ElementIndex);
        const unit = builder.addUnit(Unit.Kind.Atomic, model, SymmetryOperator.Default, elements, Unit.Trait.FastBoundary);
        units.set(label_asym_id.value(i), unit);
    }

    return units;
}

function buildCellpackAssembly(model: Model, assembly: Assembly) {
    const coordinateSystem = SymmetryOperator.create(assembly.id, Mat4.identity(), { assembly: { id: assembly.id, operId: 0, operList: [] } });
    const assembler = Structure.Builder({
        coordinateSystem,
        label: model.label,
    });

    const units = createModelChainMap(model);

    for (const g of assembly.operatorGroups) {
        for (const oper of g.operators) {
            for (const id of g.asymIds!) {
                assembler.addWithOperator(units.get(id)!, oper);
            }
        }
    }

    return assembler.getStructure();
}

type CellpackAssembly = typeof CellpackAssembly
const CellpackAssembly = PluginStateTransform.BuiltIn({
    name: 'cellpack-assembly',
    display: { name: 'Cellpack Assembly' },
    from: PSO.Molecule.Model,
    to: PSO.Molecule.Structure,
    params: {
        id: PD.Text('', { label: 'Asm Id', description: 'Assembly Id (use empty for the 1st assembly)' }),
    }
})({
    canAutoUpdate({ newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Build Structure', async ctx => {
            const model = a.data;

            let id = params.id;
            let asm: Assembly | undefined = void 0;
            const symmetry = ModelSymmetry.Provider.get(model);

            // if no id is specified, use the 1st assembly.
            if (!id && symmetry && symmetry.assemblies.length !== 0) {
                id = symmetry.assemblies[0].id;
            }

            if (!symmetry || symmetry.assemblies.length === 0) {
                plugin.log.warn(`Model '${model.entryId}' has no assembly, returning model structure.`);
            } else {
                asm = Symmetry.findAssembly(model, id || '');
                if (!asm) {
                    plugin.log.warn(`Model '${model.entryId}' has no assembly called '${id}', returning model structure.`);
                }
            }

            const base = Structure.ofModel(model);
            if (!asm) {
                const label = { label: 'Model', description: Structure.elementDescription(base) };
                return new PSO.Molecule.Structure(base, label);
            }

            const s = buildCellpackAssembly(model, asm);

            const objProps = { label: `Assembly ${id}`, description: Structure.elementDescription(s) };
            return new PSO.Molecule.Structure(s, objProps);
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

const CreateTransformer = StateTransformer.builderFactory('cellPACK');
export const CreateCompartmentSphere = CreateTransformer({
    name: 'create-compartment-sphere',
    display: 'Compartment Sphere',
    from: PSO.Root, // or whatever data source
    to: PSO.Shape.Representation3D,
    params: {
        center: PD.Vec3(Vec3()),
        radius: PD.Numeric(1),
        label: PD.Text('Compartment Sphere')
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Compartment Sphere', async ctx => {
            const data = params;
            const repr = MBRepresentation({ webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes }, () => (MBParams));
            await repr.createOrUpdate({ ...params, quality: 'custom', xrayShaded: true, doubleSided: true }, data).runInContext(ctx);
            return new PSO.Shape.Representation3D({ repr, sourceData: a }, { label: data.label });
        });
    }
});

type UnitsByEntity = Map<EntityIndex, Unit[]>;
const UnitsByEntity = CustomStructureProperty.createSimple<UnitsByEntity>('units_by_entity', 'local');

function getUnitsByEntity(structure: Structure): UnitsByEntity {
    if (UnitsByEntity.get(structure).value) {
        return UnitsByEntity.get(structure).value!;
    }

    const atomicIndex = structure.model.atomicHierarchy.index;
    const map: UnitsByEntity = new Map();
    for (const ug of structure.unitSymmetryGroups) {
        const u = ug.units[0] as Unit.Atomic;
        const e = atomicIndex.getEntityFromChain(u.chainIndex[u.elements[0]]);

        if (!map.has(e)) map.set(e, []);
        const entityUnits = map.get(e)!;

        for (let i = 0, il = ug.units.length; i < il; ++i) {
            entityUnits.push(ug.units[i]);
        }
    }

    UnitsByEntity.set(structure, { value: map });

    return map;
}

type EntityStructure = typeof EntityStructure
const EntityStructure = PluginStateTransform.BuiltIn({
    name: 'entity-structure',
    display: { name: 'Entity Structure' },
    from: PSO.Molecule.Structure,
    to: PSO.Molecule.Structure,
    params: {
        entityId: PD.Text('')
    }
})({
    canAutoUpdate({ newParams }) {
        return true;
    },
    apply({ a, params }) {
        return Task.create('Build Structure', async ctx => {
            const parent = a.data;
            const { entities } = parent.model;
            const idx = entities.getEntityIndex(params.entityId);

            const unitsByEntity = getUnitsByEntity(parent);
            const units = unitsByEntity.get(idx) || [];
            // if (!unitsByEntity.get(idx)) {
            //     console.log(entities.data.pdbx_description.value(idx));
            // }
            const structure = Structure.create(units, { parent });

            const description = entities.data.pdbx_description.value(idx)[0] || 'model';
            const label = description.split('.').at(-1) || a.label;

            return new PSO.Molecule.Structure(structure, { label, description: `${a.description}` });
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});