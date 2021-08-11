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
import { IngredientFiles, getFloatValue } from './util';
import { Asset } from '../../mol-util/assets';
import { PluginContext } from '../../mol-plugin/context';
import { CellPackInfoProvider } from './property';
import { Structure, StructureSymmetry, Unit, Model } from '../../mol-model/structure';
import { ModelSymmetry } from '../../mol-model-formats/structure/property/symmetry';
import { Vec3, Quat } from '../../mol-math/linear-algebra';
import { readFromFile } from '../../mol-util/data-source';
import { StateTransformer } from '../../mol-state';
import { MBRepresentation, MBParams } from './representation';

//export const DefaultCellPackBaseUrl = 'https://mesoscope.scripps.edu/data/cellPACK_data/cellPACK_database_1.1.0/';
export const DefaultCellPackBaseUrl = 'https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0/'
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
            modeFile: PD.File({ accept: '.bin' })
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
            const { compartments, cytoplasme } = cell;
            let iName = "";
            if(!cell.mapping_ids)cell.mapping_ids={};
            if (cytoplasme) {
                packings.push({ name: 'Cytoplasme', location: 'cytoplasme', ingredients: cytoplasme.ingredients });
                for (iName in cytoplasme.ingredients){
                    if (cytoplasme.ingredients[iName].ingtype == 'fiber') {
                        cell.mapping_ids[-(fiber_counter_id+1)]=[comp_counter,iName];
                        if (!cytoplasme.ingredients[iName].nbCurve) cytoplasme.ingredients[iName].nbCurve = 0;
                        fiber_counter_id++;                        
                    }
                    else {
                        cell.mapping_ids[counter_id]=[comp_counter,iName];
                        if (!cytoplasme.ingredients[iName].results) {cytoplasme.ingredients[iName].results=[]}
                        counter_id++;
                    }
                }
                comp_counter++;
            }
            if (compartments) {
                for (const name in compartments) {
                    const { surface, interior } = compartments[name];
                    if (surface) {
                        packings.push({ name, location: 'surface', ingredients: surface.ingredients, mb: compartments[name].mb });
                        for (iName in surface.ingredients){
                            if (surface.ingredients[iName].ingtype == 'fiber') {
                                cell.mapping_ids[-(fiber_counter_id+1)]=[comp_counter,iName];
                                if (!surface.ingredients[iName].nbCurve) surface.ingredients[iName].nbCurve = 0;
                                fiber_counter_id++;                        
                            }
                            else {
                                cell.mapping_ids[counter_id]=[comp_counter,iName];
                                if (!surface.ingredients[iName].results) {surface.ingredients[iName].results=[]}
                                counter_id++;
                            }
                        }  
                        comp_counter++;                  
                    }
                    if (interior) {
                        packings.push({ name, location: 'interior', ingredients: interior.ingredients });
                        for (iName in interior.ingredients){
                            if (interior.ingredients[iName].ingtype == 'fiber') {
                                cell.mapping_ids[-(fiber_counter_id+1)]=[comp_counter,iName];
                                if (!interior.ingredients[iName].nbCurve) interior.ingredients[iName].nbCurve = 0;
                                fiber_counter_id++;                        
                            }
                            else {
                                cell.mapping_ids[counter_id]=[comp_counter,iName];
                                if (!interior.ingredients[iName].results) {interior.ingredients[iName].results=[]}
                                counter_id++;
                            }
                        }  
                        comp_counter++;                  
                    }
                }
            }
            if (params.modeFile && params.modeFile.file){
                const model_data = await readFromFile(params.modeFile.file, 'binary').runInContext(ctx);//async ?
                var numbers = new DataView(model_data.buffer);
                var ninst = getFloatValue(numbers,0);
                var npoints = getFloatValue(numbers,4);
                var ncurve = getFloatValue(numbers,8);
                /*
                console.log("Parse CellPack");
                console.log("ninst ",ninst);
                console.log("npoints ",npoints);
                console.log("ncurve ",ncurve);
                */
                let pos = new Float32Array();
                let quat = new Float32Array();
                let ctr_pos = new Float32Array();
                //let ctr_norm = new Float32Array();
                let ctr_info = new Float32Array();
                let curve_ids = new Float32Array();                                
                //let ctrl_pts= new Float32Array();
                //let ctrl_normal= new Float32Array();
                //let ctrl_info= new Float32Array();//#curve_id, curve_type, angle, uLength
                let offset = 12;
                if (ninst !== 0){
                    pos = new Float32Array(model_data.buffer,offset,ninst*4);offset+=ninst*4*4;
                    quat = new Float32Array(model_data.buffer,offset,ninst*4);offset+=ninst*4*4;
                }
                if ( npoints != 0 )
                {
                    ctr_pos = new Float32Array(model_data.buffer,offset,npoints*4);offset+=npoints*4*4;
                    //ctr_norm = new Float32Array(model_data,offset,ncurve*4);
                    offset+=npoints*4*4;
                    ctr_info = new Float32Array(model_data.buffer,offset,npoints*4);offset+=npoints*4*4;
                    curve_ids = new Float32Array(model_data.buffer,offset,ncurve*4);offset+=ncurve*4*4;
                }
                //data_from_buffer = {"pos":pos,"quat":quat,"ctrl_pts":ctrl_pts,"ctrl_normal":ctrl_normal,"ctrl_info":ctrl_info};
                //create all the bodys then the instance ?
                for (var i=0;i<ninst;i++)
                {
                  let x:number =  pos[i*4+0];
                  let y:number =  pos[i*4+1];
                  let z:number =  pos[i*4+2];
                  //q = new THREE.Quaternion(quats[i*4+0],quats[i*4+1],quats[i*4+2],quats[i*4+3]);
                  let ingr_id = pos[i*4+3] as number;
                  let pid = cell.mapping_ids[ingr_id];
                  if (!packings[pid[0]].ingredients[pid[1]].results) {
                    packings[pid[0]].ingredients[pid[1]].results=[];
                  }
                  packings[pid[0]].ingredients[pid[1]].results.push([Vec3.create(x,y,z), 
                                                                     Quat.create(quat[i*4+0],quat[i*4+1],quat[i*4+2],quat[i*4+3])]);
                }
                let counter = 0;
                let ctr_points:Vec3[] = [];
                let prev_ctype = 0;
                let prev_cid = 0;
                /*
                console.log("ctr_info",ctr_info);
                console.log("curve_ids",curve_ids);
                */
                for (var i=0;i<npoints;i++)
                {
                    let x:number = -ctr_pos[i*4+0];
                    let y:number =  ctr_pos[i*4+1];
                    let z:number =  ctr_pos[i*4+2];
                    let cid:number = ctr_info[i*4+0];//curve id
                    let ctype:number = curve_ids[cid*4+0];//curve type  
                    //cid  148 165 -1 0
                    //console.log("cid ",cid,ctype,prev_cid,prev_ctype);//165,148
                    if (prev_ctype != ctype){
                        let pid = cell.mapping_ids[-prev_ctype-1];
                        const cname = `curve${counter}`;
                        //console.log("ctype ",pid,cname,prev_ctype,(-prev_ctype-1));
                        packings[pid[0]].ingredients[pid[1]].nbCurve = counter+1;
                        packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
                        ctr_points = [];
                        counter = 0;
                    }
                    else if (prev_cid != cid){
                        ctr_points = [];
                        let pid = cell.mapping_ids[-prev_ctype-1];
                        const cname = `curve${counter}`;
                        //console.log("cid ",pid,cname,prev_ctype,(-prev_ctype-1),prev_cid,cid);
                        packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;            
                        counter+=1;
                    }
                    ctr_points.push(Vec3.create(x,y,z));
                    prev_ctype = ctype;
                    prev_cid = cid;
                }
                //do the last one
                if ( npoints != 0 )
                {
                    let pid = cell.mapping_ids[-prev_ctype-1];
                    const cname = `curve${counter}`;
                    //console.log("ctype ",pid,cname,prev_ctype,(-prev_ctype-1));
                    packings[pid[0]].ingredients[pid[1]].nbCurve = counter+1;
                    packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
                    ctr_points = [];
                }
                counter = 0;                
                //console.log("done");                
                //console.log(packings);
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
            return new PSO.Molecule.Structure(structure, { label: packing.name+"."+packing.location });
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
                const builder = Structure.Builder({ label: "Membrane" });
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
            await repr.createOrUpdate({ ...params, quality: 'custom',doubleSided:true }, data).runInContext(ctx);
            return new PSO.Shape.Representation3D({ repr, sourceData: a }, { label: `My Sphere` });
        });
    }
});