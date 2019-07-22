/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateAction } from '../../../../mol-state';
import { PluginContext } from '../../../../mol-plugin/context';
import { PluginStateObject as PSO } from '../../../../mol-plugin/state/objects';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Ingredient, Packing } from './data';
import { getFromPdb, getFromCellPackDB } from './util';
import { Model, Structure, StructureSymmetry, StructureSelection, Queries, QueryContext, StructureProperties as SP } from '../../../../mol-model/structure';
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

async function getStructure(model: Model, props: { assembly?: string, trace?: boolean }) {
    let structure = Structure.ofModel(model)
    const { assembly, trace } = props

    if (assembly) {
        structure = await StructureSymmetry.buildAssembly(structure, assembly).run()
    }

    if (trace === true) {
        const query = Queries.generators.atoms({ atomTest: ctx => {
            const atomId = SP.atom.label_atom_id(ctx.element)
            return atomId === 'CA' || atomId === 'P'
        }})
        const result = query(new QueryContext(structure))
        structure = StructureSelection.unionStructure(result)
    }

    const query = Queries.internal.atomicSequence()
    const result = query(new QueryContext(structure))
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

function getTransforms(results: Ingredient['results']) {
    return results.map((r: Ingredient['results'][0]) => getTransform(r[0], r[1]))
}

function getAssembly(transforms: Mat4[], structure: Structure) {
    const builder = Structure.Builder(void 0, void 0)
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
    if (source.pdb === 'None') return

    // TODO handle fibers
    if (nbCurve) return

    const model = await getModel(source.pdb || name, baseUrl)
    if (!model) return

    const structure = await getStructure(model, { trace: true, assembly: source.biomt ? '1' : undefined })
    const transforms = getTransforms(results)
    const assembly = getAssembly(transforms, structure)
    return assembly
}

export function createStructureFromCellPack(ingredients: Packing['ingredients'], baseUrl: string) {
    return Task.create('Create Packing Structure', async ctx => {
        const structures: Structure[] = []
        for (const iName in ingredients) {
            if (ctx.shouldUpdate) ctx.update(iName)
            const s = await getIngredientStructure(ingredients[iName], baseUrl)
            if (s) structures.push(s)
        }

        const builder = Structure.Builder(void 0, void 0)
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
            ['BloodPlasma1.2.apr.json', 'BloodPlasma1.2'],
            ['BloodSerumfillResult.apr', 'BloodSerumfillResult'],
            ['HIV-1_0.1.6-8_mixed_radii_pdb.cpr', 'HIV-1_0.1.6-8_mixed_radii_pdb'],
            ['influenza_model1.json', 'influenza_model1'],
            ['Mycoplasma1.5_mixed_pdb_fixed.cpr', 'Mycoplasma1.5_mixed_pdb_fixed'],
            ['NM_Analysis_FigureC1.4.cpr.json', 'NM_Analysis_FigureC1.4']
        ]),
        baseUrl: PD.Text('https://cdn.jsdelivr.net/gh/mesoscope/cellPACK_data@master/')
    },
    from: PSO.Root
})(({ state, params }, ctx: PluginContext) => Task.create('CellPack Loader', async taskCtx => {
    const url = getCellPackModelUrl(params.id, params.baseUrl)

    const root = state.build().toRoot();

    const cellPackBuilder = root
        .apply(StateTransforms.Data.Download, { url, isBinary: false, label: params.id })
        .apply(StateTransforms.Data.ParseJson)
        .apply(ParseCellPack)

    const cellPackObject = await state.updateTree(cellPackBuilder).runInContext(taskCtx)
    const { packings } = cellPackObject.data
    let tree = state.build().to(cellPackBuilder.ref);

    const colors = distinctColors(packings.length)

    for (let i = 0, il = packings.length; i < il; ++i) {
        const hcl = Hcl.fromColor(Hcl(), colors[i])
        const hue = [Math.max(0, hcl[0] - 35), Math.min(360, hcl[0] + 35)] as [number, number]

        tree.apply(StructureFromCellpack, { packing: i, baseUrl: params.baseUrl })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.createParams(ctx, Structure.Empty, {
                    // repr: ctx.structureRepresentation.registry.get('point'),
                    repr: [
                        ctx.structureRepresentation.registry.get('gaussian-surface'),
                        (c, ctx) => {
                            return {
                                quality: 'custom', resolution: 10, radiusOffset: 2,
                                alpha: 1.0, flatShaded: false, doubleSided: false,
                            }
                        }
                    ],
                    color: [
                        ModelIndexColorThemeProvider,
                        (c, ctx) => {
                            return {
                                palette: {
                                    name: 'generate',
                                    params: {
                                        hue, chroma: [30, 80], luminance: [15, 85],
                                        clusteringStepCount: 50, minSampleCount: 800
                                    }
                                }
                            }
                        }
                    ]
                }))
            }

    await state.updateTree(tree).runInContext(taskCtx);
}));