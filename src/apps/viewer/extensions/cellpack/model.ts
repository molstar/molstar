/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateAction } from '../../../../mol-state';
import { PluginContext } from '../../../../mol-plugin/context';
import { PluginStateObject as PSO } from '../../../../mol-plugin/state/objects';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Ingredient, CellPacking } from './data';
import { getFromPdb, getFromCellPackDB } from './util';
import { Model, Structure, StructureSymmetry, StructureSelection, QueryContext } from '../../../../mol-model/structure';
import { trajectoryFromMmCIF } from '../../../../mol-model-formats/structure/mmcif';
import { trajectoryFromPDB } from '../../../../mol-model-formats/structure/pdb';
import { Mat4, Vec3, Quat } from '../../../../mol-math/linear-algebra';
import { SymmetryOperator } from '../../../../mol-math/geometry';
import { Task } from '../../../../mol-task';
import { StructureRepresentation3DHelpers } from '../../../../mol-plugin/state/transforms/representation';
import { StateTransforms } from '../../../../mol-plugin/state/transforms';
import { distinctColors } from '../../../../mol-util/color/distinct';
import { ModelIndexColorThemeProvider } from '../../../../mol-theme/color/model-index';
import { Hcl } from '../../../../mol-util/color/spaces/hcl';
import { ParseCellPack, StructureFromCellpack } from './state';
import { formatMolScript } from '../../../../mol-script/language/expression-formatter';
import { MolScriptBuilder as MS } from '../../../../mol-script/language/builder';
import { getMatFromResamplePoints } from './curve';
import { compile } from '../../../../mol-script/runtime/query/compiler';
import { UniformColorThemeProvider } from '../../../../mol-theme/color/uniform';
import { ThemeRegistryContext } from '../../../../mol-theme/theme';
import { ColorTheme } from '../../../../mol-theme/color';

function getCellPackModelUrl(fileName: string, baseUrl: string) {
    return `${baseUrl}/cellPACK_database_1.1.0/results/${fileName}`
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
        const op = SymmetryOperator.create(id, transforms[i], { id, operList: [ id ] })
        for (const unit of units) {
            builder.addWithOperator(unit, op)
        }
    }

    return builder.getStructure();
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

    const structure = await getStructure(model, { assembly: source.biomt ? '1' : undefined })
    const transforms = nbCurve ? getCurveTransforms(ingredient) : getResultTransforms(results)
    const assembly = getAssembly(transforms, structure)
    return assembly
}

export function createStructureFromCellPack(packing: CellPacking, baseUrl: string) {
    return Task.create('Create Packing Structure', async ctx => {
        const { ingredients, name } = packing
        const structures: Structure[] = []
        for (const iName in ingredients) {
            if (ctx.shouldUpdate) ctx.update(iName)
            const s = await getIngredientStructure(ingredients[iName], baseUrl)
            if (s) structures.push(s)
        }

        const builder = Structure.Builder({ label: name })
        let offsetInvariantId = 0
        for (const s of structures) {
            let maxInvariantId = 0
            for (const u of s.units) {
                const invariantId = u.invariantId + offsetInvariantId
                if (u.invariantId > maxInvariantId) maxInvariantId = u.invariantId
                builder.addUnit(u.kind, u.model, u.conformation.operator, u.elements, invariantId)
            }
            offsetInvariantId += maxInvariantId
        }

        const s = builder.getStructure()
        return s
    })
}

export const LoadCellPackModel = StateAction.build({
    display: { name: 'Load CellPack Model' },
    params: {
        id: PD.Select('influenza_model1.json', [
            ['blood_hiv_immature_inside.json', 'blood_hiv_immature_inside'],
            ['BloodHIV1.0_mixed_fixed_nc1.cpr', 'BloodHIV1.0_mixed_fixed_nc1'],
            ['HIV-1_0.1.6-8_mixed_radii_pdb.cpr', 'HIV-1_0.1.6-8_mixed_radii_pdb'],
            ['influenza_model1.json', 'influenza_model1'],
            ['Mycoplasma1.5_mixed_pdb_fixed.cpr', 'Mycoplasma1.5_mixed_pdb_fixed'],
        ]),
        baseUrl: PD.Text('https://cdn.jsdelivr.net/gh/mesoscope/cellPACK_data@master'),
        preset: PD.Group({
            traceOnly: PD.Boolean(false),
            representation: PD.Select('spacefill', [
                ['spacefill', 'Spacefill'],
                ['gaussian-surface', 'Gaussian Surface'],
                ['point', 'Point'],
            ])
        }, { isExpanded: true })
    },
    from: PSO.Root
})(({ state, params }, ctx: PluginContext) => Task.create('CellPack Loader', async taskCtx => {
    const url = getCellPackModelUrl(params.id, params.baseUrl)

    const root = state.build().toRoot();

    const cellPackBuilder = root
        .apply(StateTransforms.Data.Download, { url, isBinary: false, label: params.id }, { state: { isGhost: true } })
        .apply(StateTransforms.Data.ParseJson, undefined, { state: { isGhost: true } })
        .apply(ParseCellPack)

    const cellPackObject = await state.updateTree(cellPackBuilder).runInContext(taskCtx)
    const { packings } = cellPackObject.data
    let tree = state.build().to(cellPackBuilder.ref);

    const isHiv = (
        params.id === 'BloodHIV1.0_mixed_fixed_nc1.cpr' ||
        params.id === 'HIV-1_0.1.6-8_mixed_radii_pdb.cpr'
    )

    if (isHiv) {
        for (let i = 0, il = packings.length; i < il; ++i) {
            if (packings[i].name === 'HIV1_capsid_3j3q_PackInner_0_1_0') {
                const url = `${params.baseUrl}/cellPACK_database_1.1.0/extras/rna_allpoints.json`
                const data = await ctx.fetch({ url, type: 'string' }).runInContext(taskCtx);
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

    const colors = distinctColors(packings.length)

    for (let i = 0, il = packings.length; i < il; ++i) {
        const hcl = Hcl.fromColor(Hcl(), colors[i])
        const hue = [Math.max(0, hcl[0] - 35), Math.min(360, hcl[0] + 35)] as [number, number]
        const p = { packing: i, baseUrl: params.baseUrl }

        const expression = params.preset.traceOnly
            ? MS.struct.generator.atomGroups({
                'atom-test': MS.core.logic.or([
                    MS.core.rel.eq([MS.ammp('label_atom_id'), 'CA']),
                    MS.core.rel.eq([MS.ammp('label_atom_id'), 'P'])
                ])
            })
            : MS.struct.generator.all()
        const query = { language: 'mol-script' as const, expression: formatMolScript(expression) }

        tree.apply(StructureFromCellpack, p)
            .apply(StateTransforms.Model.UserStructureSelection, { query })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.createParams(ctx, Structure.Empty, {
                    repr: getReprParams(ctx, params.preset),
                    color: getColorParams(hue)
                })
            )
    }

    if (isHiv) {
        const url = `${params.baseUrl}/cellPACK_database_1.1.0/membranes/hiv_lipids.bcif`
        tree.apply(StateTransforms.Data.Download, { url }, { state: { isGhost: true } })
            .apply(StateTransforms.Data.ParseCif, undefined, { state: { isGhost: true } })
            .apply(StateTransforms.Model.TrajectoryFromMmCif, undefined, { state: { isGhost: true } })
            .apply(StateTransforms.Model.ModelFromTrajectory, undefined, { state: { isGhost: true } })
            .apply(StateTransforms.Model.StructureFromModel, undefined, { state: { isGhost: true } })
            .apply(StateTransforms.Misc.CreateGroup, { label: 'HIV1_envelope_Membrane' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.createParams(ctx, Structure.Empty, {
                    repr: getReprParams(ctx, params.preset),
                    color: UniformColorThemeProvider
                })
            )
    }

    await state.updateTree(tree).runInContext(taskCtx);
}));

function getReprParams(ctx: PluginContext, params: { representation: 'spacefill' | 'gaussian-surface' | 'point', traceOnly: boolean }) {
    const { representation, traceOnly } = params
    switch (representation) {
        case 'spacefill':
            return traceOnly
                ? [
                    ctx.structureRepresentation.registry.get('spacefill'),
                    () => ({ sizeFactor: 2, ignoreHydrogens: true })
                ] as [any, any]
                : [
                    ctx.structureRepresentation.registry.get('spacefill'),
                    () => ({ ignoreHydrogens: true })
                ] as [any, any]
        case 'gaussian-surface':
            return [
                ctx.structureRepresentation.registry.get('gaussian-surface'),
                () => ({
                    quality: 'custom', resolution: 10, radiusOffset: 2,
                    alpha: 1.0, flatShaded: false, doubleSided: false,
                    ignoreHydrogens: true
                })
            ] as [any, any]
        case 'point':
            return [
                ctx.structureRepresentation.registry.get('point'),
                () => ({ ignoreHydrogens: true })
            ] as [any, any]
    }
}

function getColorParams(hue: [number, number]) {
    return [
        ModelIndexColorThemeProvider,
        (c: ColorTheme.Provider<any>, ctx: ThemeRegistryContext) => {
            return {
                palette: {
                    name: 'generate',
                    params: {
                        hue, chroma: [30, 80], luminance: [15, 85],
                        clusteringStepCount: 50, minSampleCount: 800,
                        maxCount: 75
                    }
                }
            }
        }
    ] as [any, any]
}