/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateAction, StateBuilder, StateTransformer } from '../../../../mol-state';
import { PluginContext } from '../../../../mol-plugin/context';
import { PluginStateObject as PSO } from '../../../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Ingredient, CellPacking } from './data';
import { getFromPdb, getFromCellPackDB } from './util';
import { Model, Structure, StructureSymmetry, StructureSelection, QueryContext, Unit } from '../../../../mol-model/structure';
import { trajectoryFromMmCIF, MmcifFormat } from '../../../../mol-model-formats/structure/mmcif';
import { trajectoryFromPDB } from '../../../../mol-model-formats/structure/pdb';
import { Mat4, Vec3, Quat } from '../../../../mol-math/linear-algebra';
import { SymmetryOperator } from '../../../../mol-math/geometry';
import { Task, RuntimeContext } from '../../../../mol-task';
import { StateTransforms } from '../../../../mol-plugin-state/transforms';
import { distinctColors } from '../../../../mol-util/color/distinct';
import { Hcl } from '../../../../mol-util/color/spaces/hcl';
import { ParseCellPack, StructureFromCellpack, DefaultCellPackBaseUrl } from './state';
import { MolScriptBuilder as MS } from '../../../../mol-script/language/builder';
import { getMatFromResamplePoints } from './curve';
import { compile } from '../../../../mol-script/runtime/query/compiler';
import { CifCategory, CifField } from '../../../../mol-io/reader/cif';
import { mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif';
import { Column } from '../../../../mol-data/db';
import { createModels } from '../../../../mol-model-formats/structure/basic/parser';
import { CellpackPackingsPreset } from './preset';
import { AjaxTask } from '../../../../mol-util/data-source';

function getCellPackModelUrl(fileName: string, baseUrl: string) {
    return `${baseUrl}/results/${fileName}`
}

async function getModel(id: string, baseUrl: string) {
    let model: Model;
    if (id.match(/^[1-9][a-zA-Z0-9]{3,3}$/i)) {
        // return
        const cif = await getFromPdb(id)
        model = (await trajectoryFromMmCIF(cif).run())[0]
    } else {
        const pdb = await getFromCellPackDB(id, baseUrl)
        model = (await trajectoryFromPDB(pdb).run())[0]
    }
    return model
}

async function getStructure(model: Model, props: { assembly?: string } = {}) {
    let structure = Structure.ofModel(model)
    const { assembly } = props

    if (assembly) {
        structure = await StructureSymmetry.buildAssembly(structure, assembly).run()
    }

    const query = MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
            'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer'])
        })
    ])
    const compiled = compile<StructureSelection>(query)
    const result = compiled(new QueryContext(structure))
    structure = StructureSelection.unionStructure(result)

    return structure
}

function getTransform(trans: Vec3, rot: Quat) {
    const q: Quat = Quat.create(-rot[3], rot[0], rot[1], rot[2])
    const m: Mat4 = Mat4.fromQuat(Mat4.zero(), q)
    Mat4.transpose(m, m)
    Mat4.scale(m, m, Vec3.create(-1.0, 1.0, -1.0))
    Mat4.setTranslation(m, trans)
    return m
}

function getResultTransforms(results: Ingredient['results']) {
    return results.map((r: Ingredient['results'][0]) => getTransform(r[0], r[1]))
}

function getCurveTransforms(ingredient: Ingredient) {
    const n = ingredient.nbCurve || 0
    const instances: Mat4[] = []

    for (let i = 0; i < n; ++i) {
        const cname = `curve${i}`
        if (!(cname in ingredient)) {
            // console.warn(`Expected '${cname}' in ingredient`)
            continue
        }
        const _points = ingredient[cname] as Vec3[]
        if (_points.length <= 2) {
            // TODO handle curve with 2 or less points
            continue
        }
        const points = new Float32Array(_points.length * 3)
        for (let i = 0, il = _points.length; i < il; ++i) Vec3.toArray(_points[i], points, i * 3)
        const newInstances = getMatFromResamplePoints(points)
        instances.push(...newInstances)
    }

    return instances
}

function getAssembly(transforms: Mat4[], structure: Structure) {
    const builder = Structure.Builder()
    const { units } = structure;

    for (let i = 0, il = transforms.length; i < il; ++i) {
        const id = `${i + 1}`
        const op = SymmetryOperator.create(id, transforms[i], { assembly: { id, operId: i, operList: [ id ] } })
        for (const unit of units) {
            builder.addWithOperator(unit, op)
        }
    }

    return builder.getStructure();
}

function getCifCurve(name: string, transforms: Mat4[], model: Model) {
    if (!MmcifFormat.is(model.sourceData)) throw new Error('mmcif source data needed')

    const { db } = model.sourceData.data
    const d = db.atom_site
    const n = d._rowCount
    const rowCount = n * transforms.length

    const { offsets, count } = model.atomicHierarchy.chainAtomSegments

    const x = d.Cartn_x.toArray()
    const y = d.Cartn_y.toArray()
    const z = d.Cartn_z.toArray()

    const Cartn_x = new Float32Array(rowCount)
    const Cartn_y = new Float32Array(rowCount)
    const Cartn_z = new Float32Array(rowCount)
    const map = new Uint32Array(rowCount)
    const seq = new Int32Array(rowCount)
    let offset = 0
    for (let c = 0; c < count; ++c) {
        const cStart = offsets[c]
        const cEnd = offsets[c + 1]
        const cLength = cEnd - cStart
        for (let t = 0, tl = transforms.length; t < tl; ++t) {
            const m = transforms[t]
            for (let j = cStart; j < cEnd; ++j) {
                const i = offset + j - cStart
                const xj = x[j], yj = y[j], zj = z[j]
                Cartn_x[i] = m[0] * xj + m[4] * yj + m[8] * zj + m[12]
                Cartn_y[i] = m[1] * xj + m[5] * yj + m[9] * zj + m[13]
                Cartn_z[i] = m[2] * xj + m[6] * yj + m[10] * zj + m[14]
                map[i] = j
                seq[i] = t + 1
            }
            offset += cLength
        }
    }

    function multColumn<T>(column: Column<T>) {
        const array = column.toArray()
        return Column.ofLambda({
            value: row => array[map[row]],
            areValuesEqual: (rowA, rowB) => map[rowA] === map[rowB] || array[map[rowA]] === array[map[rowB]],
            rowCount, schema: column.schema
        })
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
    }

    const categories = {
        entity: CifCategory.ofTable('entity', db.entity),
        chem_comp: CifCategory.ofTable('chem_comp', db.chem_comp),
        atom_site: CifCategory.ofFields('atom_site', _atom_site)
    }

    return {
        header: name,
        categoryNames: Object.keys(categories),
        categories
    };
}

async function getCurve(name: string, transforms: Mat4[], model: Model) {
    const cif = getCifCurve(name, transforms, model)

    const curveModelTask = Task.create('Curve Model', async ctx => {
        const format = MmcifFormat.fromFrame(cif)
        const models = await createModels(format.data.db, format, ctx)
        return models[0]
    })

    const curveModel = await curveModelTask.run()
    return getStructure(curveModel)
}

async function getIngredientStructure(ingredient: Ingredient, baseUrl: string) {
    const { name, source, results, nbCurve } = ingredient

    // TODO can these be added to the library?
    if (name === 'HIV1_CAhex_0_1_0') return
    if (name === 'HIV1_CAhexCyclophilA_0_1_0') return
    if (name === 'iLDL') return
    if (name === 'peptides') return
    if (name === 'lypoglycane') return

    if (source.pdb === 'None') return

    const model = await getModel(source.pdb || name, baseUrl)
    if (!model) return

    if (nbCurve) {
        return getCurve(name, getCurveTransforms(ingredient), model)
    } else {
        const structure = await getStructure(model, { assembly: source.biomt ? '1' : undefined })
        return getAssembly(getResultTransforms(results), structure)
    }
}

export function createStructureFromCellPack(packing: CellPacking, baseUrl: string) {
    return Task.create('Create Packing Structure', async ctx => {
        const { ingredients, name } = packing
        const structures: Structure[] = []
        for (const iName in ingredients) {
            if (ctx.shouldUpdate) await ctx.update(iName)
            const s = await getIngredientStructure(ingredients[iName], baseUrl)
            if (s) structures.push(s)
        }

        if (ctx.shouldUpdate) await ctx.update(`${name} - units`)
        const builder = Structure.Builder({ label: name })
        let offsetInvariantId = 0
        for (const s of structures) {
            if (ctx.shouldUpdate) await ctx.update(`${s.label}`)
            let maxInvariantId = 0
            for (const u of s.units) {
                const invariantId = u.invariantId + offsetInvariantId
                if (u.invariantId > maxInvariantId) maxInvariantId = u.invariantId
                builder.addUnit(u.kind, u.model, u.conformation.operator, u.elements, Unit.Trait.None, invariantId)
            }
            offsetInvariantId += maxInvariantId + 1
        }

        if (ctx.shouldUpdate) await ctx.update(`${name} - structure`)
        const s = builder.getStructure()
        for( let i = 0, il = s.models.length; i < il; ++i) {
            const { trajectoryInfo } = s.models[i]
            trajectoryInfo.size = il
            trajectoryInfo.index = i
        }
        return s
    })
}

async function handleHivRna(ctx: { runtime: RuntimeContext, fetch: AjaxTask }, packings: CellPacking[], baseUrl: string) {
    for (let i = 0, il = packings.length; i < il; ++i) {
        if (packings[i].name === 'HIV1_capsid_3j3q_PackInner_0_1_0') {
            const url = `${baseUrl}/extras/rna_allpoints.json`
            const data = await ctx.fetch({ url, type: 'string' }).runInContext(ctx.runtime);
            const { points } = await (new Response(data)).json() as { points: number[] }

            const curve0: Vec3[] = []
            for (let j = 0, jl = points.length; j < jl; j += 3) {
                curve0.push(Vec3.fromArray(Vec3(), points, j))
            }
            packings[i].ingredients['RNA'] = {
                source: { pdb: 'RNA_U_Base.pdb', transform: { center: false } },
                results: [],
                name: 'RNA',
                nbCurve: 1,
                curve0
            }
        }
    }
}

export const LoadCellPackModel = StateAction.build({
    display: { name: 'Load CellPack Model' },
    params: {
        source: PD.MappedStatic('id', {
            'id': PD.Select('influenza_model1.json', [
                ['blood_hiv_immature_inside.json', 'blood_hiv_immature_inside'],
                ['BloodHIV1.0_mixed_fixed_nc1.cpr', 'BloodHIV1.0_mixed_fixed_nc1'],
                ['HIV-1_0.1.6-8_mixed_radii_pdb.cpr', 'HIV-1_0.1.6-8_mixed_radii_pdb'],
                ['influenza_model1.json', 'influenza_model1'],
                ['Mycoplasma1.5_mixed_pdb_fixed.cpr', 'Mycoplasma1.5_mixed_pdb_fixed'],
            ] as const),
            'file': PD.File({ accept: 'id' }),
        }, { options: [['id', 'Id'], ['file', 'File']] }),
        // id: PD.Select('influenza_model1.json', [
        //     ['blood_hiv_immature_inside.json', 'blood_hiv_immature_inside'],
        //     ['BloodHIV1.0_mixed_fixed_nc1.cpr', 'BloodHIV1.0_mixed_fixed_nc1'],
        //     ['HIV-1_0.1.6-8_mixed_radii_pdb.cpr', 'HIV-1_0.1.6-8_mixed_radii_pdb'],
        //     ['influenza_model1.json', 'influenza_model1'],
        //     ['Mycoplasma1.5_mixed_pdb_fixed.cpr', 'Mycoplasma1.5_mixed_pdb_fixed'],
        //     ['curveTest', 'Curve Test'],
        // ] as const),
        baseUrl: PD.Text(DefaultCellPackBaseUrl),
        preset: PD.Group({
            traceOnly: PD.Boolean(false),
            representation: PD.Select('spacefill', PD.arrayToOptions(['spacefill', 'gaussian-surface', 'point', 'ellipsoid']))
        }, { isExpanded: true })
    },
    from: PSO.Root
})(({ state, params }, ctx: PluginContext) => Task.create('CellPack Loader', async taskCtx => {
    let cellPackJson: StateBuilder.To<PSO.Format.Json, StateTransformer<PSO.Data.String, PSO.Format.Json>>
    if (params.source.name === 'id') {
        const url = getCellPackModelUrl(params.source.params, params.baseUrl)
        cellPackJson = state.build().toRoot()
            .apply(StateTransforms.Data.Download, { url, isBinary: false, label: params.source.params }, { state: { isGhost: true } })
    } else {
        const file = params.source.params
        cellPackJson = state.build().toRoot()
            .apply(StateTransforms.Data.ReadFile, { file, isBinary: false, label: file.name }, { state: { isGhost: true } })
    }

    const cellPackBuilder = cellPackJson
        .apply(StateTransforms.Data.ParseJson, undefined, { state: { isGhost: true } })
        .apply(ParseCellPack)

    const cellPackObject = await state.updateTree(cellPackBuilder).runInContext(taskCtx)
    const { packings } = cellPackObject.data

    // TODO make configurable
    // const isHiv = params.source.name === 'id' && (
    //     params.source.params === 'BloodHIV1.0_mixed_fixed_nc1.cpr' ||
    //     params.source.params === 'HIV-1_0.1.6-8_mixed_radii_pdb.cpr'
    // )
    const isHiv = false

    if (isHiv) {
        await handleHivRna({ runtime: taskCtx, fetch: ctx.fetch }, packings, params.baseUrl)
    }

    const colors = distinctColors(packings.length)
    for (let i = 0, il = packings.length; i < il; ++i) {
        const hcl = Hcl.fromColor(Hcl(), colors[i])
        const hue = [Math.max(0, hcl[0] - 35), Math.min(360, hcl[0] + 35)] as [number, number]
        const p = { packing: i, baseUrl: params.baseUrl }

        const packing = state.build().to(cellPackBuilder.ref).apply(StructureFromCellpack, p)
        await ctx.updateDataState(packing, { revertOnError: true });

        const packingParams = {
            traceOnly: params.preset.traceOnly,
            representation: params.preset.representation,
            hue
        }
        await CellpackPackingsPreset.apply(packing.selector, packingParams, ctx)
    }

    // TODO
    // if (isHiv) {
    //     const url = `${params.baseUrl}/membranes/hiv_lipids.bcif`
    //     tree.apply(StateTransforms.Data.Download, { label: 'hiv_lipids', url, isBinary: true }, { state: { isGhost: true } })
    //         .apply(StateTransforms.Data.ParseCif, undefined, { state: { isGhost: true } })
    //         .apply(StateTransforms.Model.TrajectoryFromMmCif, undefined, { state: { isGhost: true } })
    //         .apply(StateTransforms.Model.ModelFromTrajectory, undefined, { state: { isGhost: true } })
    //         .apply(StateTransforms.Model.StructureFromModel, undefined, { state: { isGhost: true } })
    //         .apply(StateTransforms.Misc.CreateGroup, { label: 'HIV1_envelope_Membrane' })
    //         .apply(StateTransforms.Representation.StructureRepresentation3D,
    //             createStructureRepresentationParams(ctx, Structure.Empty, {
    //                 ...getReprParams(ctx, params.preset),
    //                 color: UniformColorThemeProvider
    //             })
    //         )
    // }
}));