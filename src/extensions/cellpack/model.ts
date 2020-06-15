/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateAction, StateBuilder, StateTransformer, State } from '../../mol-state';
import { PluginContext } from '../../mol-plugin/context';
import { PluginStateObject as PSO } from '../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Ingredient, IngredientSource, CellPacking } from './data';
import { getFromPdb, getFromCellPackDB, IngredientFiles, parseCif, parsePDBfile, getStructureMean, getFromOPM } from './util';
import { Model, Structure, StructureSymmetry, StructureSelection, QueryContext, Unit } from '../../mol-model/structure';
import { trajectoryFromMmCIF, MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { trajectoryFromPDB } from '../../mol-model-formats/structure/pdb';
import { Mat4, Vec3, Quat } from '../../mol-math/linear-algebra';
import { SymmetryOperator } from '../../mol-math/geometry';
import { Task, RuntimeContext } from '../../mol-task';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { ParseCellPack, StructureFromCellpack, DefaultCellPackBaseUrl, StructureFromAssemblies } from './state';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { getMatFromResamplePoints } from './curve';
import { compile } from '../../mol-script/runtime/query/compiler';
import { CifCategory, CifField } from '../../mol-io/reader/cif';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { Column } from '../../mol-data/db';
import { createModels } from '../../mol-model-formats/structure/basic/parser';
import { CellpackPackingPreset, CellpackMembranePreset } from './preset';
import { Asset } from '../../mol-util/assets';
import { Color } from '../../mol-util/color';
import { readFromFile } from '../../mol-util/data-source';
import { objectForEach } from '../../mol-util/object';

function getCellPackModelUrl(fileName: string, baseUrl: string) {
    return `${baseUrl}/results/${fileName}`;
}

class TrajectoryCache {
    private map = new Map<string, Model.Trajectory>();
    set(id: string, trajectory: Model.Trajectory) { this.map.set(id, trajectory); }
    get(id: string) { return this.map.get(id); }
}

async function getModel(plugin: PluginContext, id: string, ingredient: Ingredient, baseUrl: string, trajCache: TrajectoryCache, file?: Asset.File) {
    const assetManager = plugin.managers.asset;
    const modelIndex = (ingredient.source.model) ? parseInt(ingredient.source.model) : 0;
    const surface = (ingredient.ingtype) ? (ingredient.ingtype === 'transmembrane') : false;
    let trajectory = trajCache.get(id);
    let assets: Asset.Wrapper[] = [];
    if (!trajectory) {
        if (file) {
            if (file.name.endsWith('.cif')) {
                const text = await plugin.runTask(assetManager.resolve(file, 'string'));
                assets.push(text);
                const cif = (await parseCif(plugin, text.data)).blocks[0];
                trajectory = await plugin.runTask(trajectoryFromMmCIF(cif));
            } else if (file.name.endsWith('.bcif')) {
                const binary = await plugin.runTask(assetManager.resolve(file, 'binary'));
                assets.push(binary);
                const cif = (await parseCif(plugin, binary.data)).blocks[0];
                trajectory = await plugin.runTask(trajectoryFromMmCIF(cif));
            } else if (file.name.endsWith('.pdb')) {
                const text = await plugin.runTask(assetManager.resolve(file, 'string'));
                assets.push(text);
                const pdb = await parsePDBfile(plugin, text.data, id);
                trajectory = await plugin.runTask(trajectoryFromPDB(pdb));
            } else {
                throw new Error(`unsupported file type '${file.name}'`);
            }
        } else if (id.match(/^[1-9][a-zA-Z0-9]{3,3}$/i)) {
            if (surface){
                try {
                    const data = await getFromOPM(plugin, id, assetManager);
                    assets.push(data.asset);
                    trajectory = await plugin.runTask(trajectoryFromPDB(data.pdb));
                } catch (e) {
                    // fallback to getFromPdb
                    // console.error(e);
                    const { mmcif, asset } = await getFromPdb(plugin, id, assetManager);
                    assets.push(asset);
                    trajectory = await plugin.runTask(trajectoryFromMmCIF(mmcif));
                }
            } else {
                const { mmcif, asset } = await getFromPdb(plugin, id, assetManager);
                assets.push(asset);
                trajectory = await plugin.runTask(trajectoryFromMmCIF(mmcif));
            }
        } else {
            const data = await getFromCellPackDB(plugin, id, baseUrl, assetManager);
            assets.push(data.asset);
            if ('pdb' in data) {
                trajectory = await plugin.runTask(trajectoryFromPDB(data.pdb));
            } else {
                trajectory = await plugin.runTask(trajectoryFromMmCIF(data.mmcif));
            }
        }
        trajCache.set(id, trajectory);
    }
    const model = trajectory[modelIndex];
    return { model, assets };
}

async function getStructure(plugin: PluginContext, model: Model, source: IngredientSource, props: { assembly?: string } = {}) {
    let structure = Structure.ofModel(model);
    const { assembly } = props;

    if (assembly) {
        structure = await plugin.runTask(StructureSymmetry.buildAssembly(structure, assembly));
    }
    let query;
    if (source.selection){
        const asymIds: string[] = source.selection.replace(' ', '').replace(':', '').split('or');
        query = MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
                'chain-test': MS.core.set.has([MS.set(...asymIds), MS.ammp('auth_asym_id')])
            })
        ]);
    } else {
        query = MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer'])
            })
        ]);
    }

    const compiled = compile<StructureSelection>(query);
    const result = compiled(new QueryContext(structure));
    structure = StructureSelection.unionStructure(result);

    return structure;
}

function getTransformLegacy(trans: Vec3, rot: Quat) {
    const q: Quat = Quat.create(-rot[3], rot[0], rot[1], rot[2]);
    const m: Mat4 = Mat4.fromQuat(Mat4.zero(), q);
    Mat4.transpose(m, m);
    Mat4.scale(m, m, Vec3.create(-1.0, 1.0, -1.0));
    Mat4.setTranslation(m, trans);
    return m;
}

function getTransform(trans: Vec3, rot: Quat) {
    const q: Quat = Quat.create(rot[0], rot[1], rot[2], rot[3]);
    const m: Mat4 = Mat4.fromQuat(Mat4.zero(), q);
    const p: Vec3 = Vec3.create(trans[0], trans[1], trans[2]);
    Mat4.setTranslation(m, p);
    return m;
}

function getResultTransforms(results: Ingredient['results'], legacy: boolean) {
    if (legacy) return results.map((r: Ingredient['results'][0]) => getTransformLegacy(r[0], r[1]));
    else return results.map((r: Ingredient['results'][0]) => getTransform(r[0], r[1]));
}

function getCurveTransforms(ingredient: Ingredient) {
    const n = ingredient.nbCurve || 0;
    const instances: Mat4[] = [];
    let segmentLength = 3.4;
    if (ingredient.uLength){
        segmentLength = ingredient.uLength;
    } else if (ingredient.radii){
        segmentLength = ingredient.radii[0].radii
            ? ingredient.radii[0].radii[0] * 2.0
            : 3.4;
    }
    let resampling: boolean = false;
    for (let i = 0; i < n; ++i) {
        const cname = `curve${i}`;
        if (!(cname in ingredient)) {
            // console.warn(`Expected '${cname}' in ingredient`)
            continue;
        }
        const _points = ingredient[cname] as Vec3[];
        if (_points.length <= 2) {
            // TODO handle curve with 2 or less points
            continue;
        }
        // test for resampling
        let distance: number = Vec3.distance(_points[0], _points[1]);
        if (distance >= segmentLength + 2.0) {
            console.info(distance);
            resampling = true;
        }
        const points = new Float32Array(_points.length * 3);
        for (let i = 0, il = _points.length; i < il; ++i) Vec3.toArray(_points[i], points, i * 3);
        const newInstances = getMatFromResamplePoints(points, segmentLength, resampling);
        instances.push(...newInstances);
    }
    return instances;
}

function getAssembly(transforms: Mat4[], structure: Structure) {
    const builder = Structure.Builder();
    const { units } = structure;

    for (let i = 0, il = transforms.length; i < il; ++i) {
        const id = `${i + 1}`;
        const op = SymmetryOperator.create(id, transforms[i], { assembly: { id, operId: i, operList: [ id ] } });
        for (const unit of units) {
            builder.addWithOperator(unit, op);
        }
    }

    return builder.getStructure();
}

function getCifCurve(name: string, transforms: Mat4[], model: Model) {
    if (!MmcifFormat.is(model.sourceData)) throw new Error('mmcif source data needed');

    const { db } = model.sourceData.data;
    const d = db.atom_site;
    const n = d._rowCount;
    const rowCount = n * transforms.length;

    const { offsets, count } = model.atomicHierarchy.chainAtomSegments;

    const x = d.Cartn_x.toArray();
    const y = d.Cartn_y.toArray();
    const z = d.Cartn_z.toArray();

    const Cartn_x = new Float32Array(rowCount);
    const Cartn_y = new Float32Array(rowCount);
    const Cartn_z = new Float32Array(rowCount);
    const map = new Uint32Array(rowCount);
    const seq = new Int32Array(rowCount);
    let offset = 0;
    for (let c = 0; c < count; ++c) {
        const cStart = offsets[c];
        const cEnd = offsets[c + 1];
        const cLength = cEnd - cStart;
        for (let t = 0, tl = transforms.length; t < tl; ++t) {
            const m = transforms[t];
            for (let j = cStart; j < cEnd; ++j) {
                const i = offset + j - cStart;
                const xj = x[j], yj = y[j], zj = z[j];
                Cartn_x[i] = m[0] * xj + m[4] * yj + m[8] * zj + m[12];
                Cartn_y[i] = m[1] * xj + m[5] * yj + m[9] * zj + m[13];
                Cartn_z[i] = m[2] * xj + m[6] * yj + m[10] * zj + m[14];
                map[i] = j;
                seq[i] = t + 1;
            }
            offset += cLength;
        }
    }

    function multColumn<T>(column: Column<T>) {
        const array = column.toArray();
        return Column.ofLambda({
            value: row => array[map[row]],
            areValuesEqual: (rowA, rowB) => map[rowA] === map[rowB] || array[map[rowA]] === array[map[rowB]],
            rowCount, schema: column.schema
        });
    }

    const _atom_site: CifCategory.SomeFields<mmCIF_Schema['atom_site']> = {
        auth_asym_id: CifField.ofColumn(multColumn(d.auth_asym_id)),
        auth_atom_id: CifField.ofColumn(multColumn(d.auth_atom_id)),
        auth_comp_id: CifField.ofColumn(multColumn(d.auth_comp_id)),
        auth_seq_id: CifField.ofNumbers(seq),

        B_iso_or_equiv: CifField.ofColumn(Column.ofConst(0, rowCount, Column.Schema.float)),
        Cartn_x: CifField.ofNumbers(Cartn_x),
        Cartn_y: CifField.ofNumbers(Cartn_y),
        Cartn_z: CifField.ofNumbers(Cartn_z),
        group_PDB: CifField.ofColumn(Column.ofConst('ATOM', rowCount, Column.Schema.str)),
        id: CifField.ofColumn(Column.ofLambda({
            value: row => row,
            areValuesEqual: (rowA, rowB) => rowA === rowB,
            rowCount, schema: d.id.schema,
        })),

        label_alt_id: CifField.ofColumn(multColumn(d.label_alt_id)),

        label_asym_id: CifField.ofColumn(multColumn(d.label_asym_id)),
        label_atom_id: CifField.ofColumn(multColumn(d.label_atom_id)),
        label_comp_id: CifField.ofColumn(multColumn(d.label_comp_id)),
        label_seq_id: CifField.ofNumbers(seq),
        label_entity_id: CifField.ofColumn(Column.ofConst('1', rowCount, Column.Schema.str)),

        occupancy: CifField.ofColumn(Column.ofConst(1, rowCount, Column.Schema.float)),
        type_symbol: CifField.ofColumn(multColumn(d.type_symbol)),

        pdbx_PDB_ins_code: CifField.ofColumn(Column.ofConst('', rowCount, Column.Schema.str)),
        pdbx_PDB_model_num: CifField.ofColumn(Column.ofConst(1, rowCount, Column.Schema.int)),
    };

    const categories = {
        entity: CifCategory.ofTable('entity', db.entity),
        chem_comp: CifCategory.ofTable('chem_comp', db.chem_comp),
        atom_site: CifCategory.ofFields('atom_site', _atom_site)
    };

    return {
        header: name,
        categoryNames: Object.keys(categories),
        categories
    };
}

async function getCurve(plugin: PluginContext, name: string, ingredient: Ingredient, transforms: Mat4[], model: Model) {
    const cif = getCifCurve(name, transforms, model);
    const curveModelTask = Task.create('Curve Model', async ctx => {
        const format = MmcifFormat.fromFrame(cif);
        const models = await createModels(format.data.db, format, ctx);
        return models[0];
    });

    const curveModel = await plugin.runTask(curveModelTask);
    return getStructure(plugin, curveModel, ingredient.source);
}

async function getIngredientStructure(plugin: PluginContext, ingredient: Ingredient, baseUrl: string, ingredientFiles: IngredientFiles, trajCache: TrajectoryCache) {
    const { name, source, results, nbCurve } = ingredient;
    if (source.pdb === 'None') return;

    const file = ingredientFiles[source.pdb];
    if (!file) {
        // TODO can these be added to the library?
        if (name === 'HIV1_CAhex_0_1_0') return; // 1VU4CtoH_hex.pdb
        if (name === 'HIV1_CAhexCyclophilA_0_1_0') return; // 1AK4fitTo1VU4hex.pdb
        if (name === 'iLDL') return; // EMD-5239
        if (name === 'peptides') return; // peptide.pdb
        if (name === 'lypoglycane') return;
    }

    // model id in case structure is NMR
    const { model, assets } = await getModel(plugin, source.pdb || name, ingredient, baseUrl, trajCache, file);
    if (!model) return;

    let structure: Structure;
    if (nbCurve) {
        structure = await getCurve(plugin, name, ingredient, getCurveTransforms(ingredient), model);
    } else {
        let bu: string|undefined = source.bu ? source.bu : undefined;
        if (bu){
            if (bu === 'AU') {
                bu = undefined;
            } else {
                bu = bu.slice(2);
            }
        }
        structure = await getStructure(plugin, model, source, { assembly: bu });
        // transform with offset and pcp
        let legacy: boolean = true;
        if (ingredient.offset || ingredient.principalAxis){
            legacy = false;
            const structureMean = getStructureMean(structure);
            Vec3.negate(structureMean, structureMean);
            const m1: Mat4 = Mat4.identity();
            Mat4.setTranslation(m1, structureMean);
            structure = Structure.transform(structure, m1);
            if (ingredient.offset){
                if (!Vec3.exactEquals(ingredient.offset, Vec3.zero())){
                    const m: Mat4 = Mat4.identity();
                    Mat4.setTranslation(m, ingredient.offset);
                    structure = Structure.transform(structure, m);
                }
            }
            if (ingredient.principalAxis){
                if (!Vec3.exactEquals(ingredient.principalAxis, Vec3.unitZ)){
                    const q: Quat = Quat.identity();
                    Quat.rotationTo(q, ingredient.principalAxis, Vec3.unitZ);
                    const m: Mat4 = Mat4.fromQuat(Mat4.zero(), q);
                    structure = Structure.transform(structure, m);
                }
            }
        }
        structure = getAssembly(getResultTransforms(results, legacy), structure);
    }

    return { structure, assets };
}

export function createStructureFromCellPack(plugin: PluginContext, packing: CellPacking, baseUrl: string, ingredientFiles: IngredientFiles) {
    return Task.create('Create Packing Structure', async ctx => {
        const { ingredients, name } = packing;
        const assets: Asset.Wrapper[] = [];
        const trajCache = new TrajectoryCache();
        const structures: Structure[] = [];
        const colors: Color[] = [];
        let skipColors: boolean = false;
        for (const iName in ingredients) {
            if (ctx.shouldUpdate) await ctx.update(iName);
            const ingredientStructure = await getIngredientStructure(plugin, ingredients[iName], baseUrl, ingredientFiles, trajCache);
            if (ingredientStructure) {
                structures.push(ingredientStructure.structure);
                assets.push(...ingredientStructure.assets);
                const c = ingredients[iName].color;
                if (c){
                    colors.push(Color.fromNormalizedRgb(c[0], c[1], c[2]));
                } else {
                    skipColors = true;
                }
            }
        }

        if (ctx.shouldUpdate) await ctx.update(`${name} - units`);
        const builder = Structure.Builder({ label: name });
        let offsetInvariantId = 0;
        for (const s of structures) {
            if (ctx.shouldUpdate) await ctx.update(`${s.label}`);
            let maxInvariantId = 0;
            for (const u of s.units) {
                const invariantId = u.invariantId + offsetInvariantId;
                if (u.invariantId > maxInvariantId) maxInvariantId = u.invariantId;
                builder.addUnit(u.kind, u.model, u.conformation.operator, u.elements, Unit.Trait.None, invariantId);
            }
            offsetInvariantId += maxInvariantId + 1;
        }

        if (ctx.shouldUpdate) await ctx.update(`${name} - structure`);
        const structure = builder.getStructure();
        for( let i = 0, il = structure.models.length; i < il; ++i) {
            Model.TrajectoryInfo.set(structure.models[i], { size: il, index: i });
        }
        return { structure, assets, colors: skipColors ? undefined : colors };
    });
}

async function handleHivRna(plugin: PluginContext, packings: CellPacking[], baseUrl: string) {
    for (let i = 0, il = packings.length; i < il; ++i) {
        if (packings[i].name === 'HIV1_capsid_3j3q_PackInner_0_1_0') {
            const url = Asset.getUrlAsset(plugin.managers.asset, `${baseUrl}/extras/rna_allpoints.json`);
            const json = await plugin.runTask(plugin.managers.asset.resolve(url, 'json', false));
            const points = json.data.points as number[];

            const curve0: Vec3[] = [];
            for (let j = 0, jl = points.length; j < jl; j += 3) {
                curve0.push(Vec3.fromArray(Vec3(), points, j));
            }
            packings[i].ingredients['RNA'] = {
                source: { pdb: 'RNA_U_Base.pdb', transform: { center: false } },
                results: [],
                name: 'RNA',
                nbCurve: 1,
                curve0
            };
        }
    }
}

async function loadMembrane(plugin: PluginContext, name: string, state: State, params: LoadCellPackModelParams) {
    let file: Asset.File | undefined = undefined;
    if (params.ingredients.files !== null) {
        const fileName = `${name}.bcif`;
        for (const f of params.ingredients.files) {
            if (fileName === f.name) {
                file = f;
                break;
            }
        }
        if (!file){
            // check for cif directly
            const cifileName = `${name}.cif`;
            for (const f of params.ingredients.files) {
                if (cifileName === f.name) {
                    file = f;
                    break;
                }
            }
        }
    }

    let b = state.build().toRoot();
    if (file) {
        if (file.name.endsWith('.cif')) {
            b = b.apply(StateTransforms.Data.ReadFile, { file, isBinary: false, label: file.name }, { state: { isGhost: true } });
        } else if (file.name.endsWith('.bcif')) {
            b = b.apply(StateTransforms.Data.ReadFile, { file, isBinary: true, label: file.name }, { state: { isGhost: true } });
        }
    } else {
        const url = Asset.getUrlAsset(plugin.managers.asset, `${params.baseUrl}/membranes/${name}.bcif`);
        b = b.apply(StateTransforms.Data.Download, { url, isBinary: true, label: name }, { state: { isGhost: true } });
    }

    const membrane = await b.apply(StateTransforms.Data.ParseCif, undefined, { state: { isGhost: true } })
        .apply(StateTransforms.Model.TrajectoryFromMmCif, undefined, { state: { isGhost: true } })
        .apply(StateTransforms.Model.ModelFromTrajectory, undefined, { state: { isGhost: true } })
        .apply(StructureFromAssemblies, undefined, { state: { isGhost: true } })
        .commit({ revertOnError: true });

    const membraneParams = {
        representation: params.preset.representation,
    };

    await CellpackMembranePreset.apply(membrane, membraneParams, plugin);
}

async function loadPackings(plugin: PluginContext, runtime: RuntimeContext, state: State, params: LoadCellPackModelParams) {
    const ingredientFiles = params.ingredients.files || [];

    let cellPackJson: StateBuilder.To<PSO.Format.Json, StateTransformer<PSO.Data.String, PSO.Format.Json>>;
    if (params.source.name === 'id') {
        const url = Asset.getUrlAsset(plugin.managers.asset, getCellPackModelUrl(params.source.params, params.baseUrl));
        cellPackJson = state.build().toRoot()
            .apply(StateTransforms.Data.Download, { url, isBinary: false, label: params.source.params }, { state: { isGhost: true } });
    } else {
        const file = params.source.params;
        if (!file?.file) {
            plugin.log.error('No file selected');
            return;
        }

        let jsonFile: Asset.File;
        if (file.name.toLowerCase().endsWith('.zip')) {
            const data = await readFromFile(file.file, 'zip').runInContext(runtime);
            jsonFile = Asset.File(new File([data['model.json']], 'model.json'));
            objectForEach(data, (v, k) => {
                if (k === 'model.json') return;
                ingredientFiles.push(Asset.File(new File([v], k)));
            });
        } else {
            jsonFile = file;
        }

        cellPackJson = state.build().toRoot()
            .apply(StateTransforms.Data.ReadFile, { file: jsonFile, isBinary: false, label: jsonFile.name }, { state: { isGhost: true } });
    }

    const cellPackBuilder = cellPackJson
        .apply(StateTransforms.Data.ParseJson, undefined, { state: { isGhost: true } })
        .apply(ParseCellPack);

    const cellPackObject = await state.updateTree(cellPackBuilder).runInContext(runtime);
    const { packings } = cellPackObject.obj!.data;

    await handleHivRna(plugin, packings, params.baseUrl);

    for (let i = 0, il = packings.length; i < il; ++i) {
        const p = { packing: i, baseUrl: params.baseUrl, ingredientFiles };

        const packing = await state.build()
            .to(cellPackBuilder.ref)
            .apply(StructureFromCellpack, p)
            .commit({ revertOnError: true });

        const packingParams = {
            traceOnly: params.preset.traceOnly,
            representation: params.preset.representation,
        };
        await CellpackPackingPreset.apply(packing, packingParams, plugin);
        if ( packings[i].location === 'surface' && params.membrane){
            await loadMembrane(plugin, packings[i].name, state, params);
        }
    }
}

const LoadCellPackModelParams = {
    source: PD.MappedStatic('id', {
        'id': PD.Select('InfluenzaModel2.json', [
            ['blood_hiv_immature_inside.json', 'Blood HIV immature'],
            ['HIV_immature_model.json', 'HIV immature'],
            ['BloodHIV1.0_mixed_fixed_nc1.cpr', 'Blood HIV'],
            ['HIV-1_0.1.6-8_mixed_radii_pdb.cpr', 'HIV'],
            ['influenza_model1.json', 'Influenza envelope'],
            ['InfluenzaModel2.json', 'Influenza Complete'],
            ['ExosomeModel.json', 'Exosome Model'],
            ['Mycoplasma1.5_mixed_pdb_fixed.cpr', 'Mycoplasma simple'],
            ['MycoplasmaModel.json', 'Mycoplasma WholeCell model'],
        ] as const, { description: 'Download the model definition with `id` from the server at `baseUrl.`' }),
        'file': PD.File({ accept: '.json,.cpr,.zip', description: 'Open model definition from .json/.cpr file or open .zip file containing model definition plus ingredients.' }),
    }, { options: [['id', 'Id'], ['file', 'File']] }),
    baseUrl: PD.Text(DefaultCellPackBaseUrl),
    membrane: PD.Boolean(true),
    ingredients : PD.Group({
        files: PD.FileList({ accept: '.cif,.bcif,.pdb' })
    }, { isExpanded: true }),
    preset: PD.Group({
        traceOnly: PD.Boolean(false),
        representation: PD.Select('gaussian-surface', PD.arrayToOptions(['spacefill', 'gaussian-surface', 'point', 'orientation']))
    }, { isExpanded: true })
};
type LoadCellPackModelParams = PD.Values<typeof LoadCellPackModelParams>

export const LoadCellPackModel = StateAction.build({
    display: { name: 'Load CellPack', description: 'Open or download a model' },
    params: LoadCellPackModelParams,
    from: PSO.Root
})(({ state, params }, ctx: PluginContext) => Task.create('CellPack Loader', async taskCtx => {
    await loadPackings(ctx, taskCtx, state, params);
}));