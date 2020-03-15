/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO: move to an "example"

import { PluginContext } from '../../mol-plugin/context';
import { Mat4 } from '../../mol-math/linear-algebra';
import { StateHelper } from './helpers';
import { PluginCommands } from '../../mol-plugin/commands';
import { StateSelection, StateBuilder } from '../../mol-state';
import { PluginStateObject as PSO } from '../../mol-plugin-state/objects';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { compile } from '../../mol-script/runtime/query/compiler';
import { StructureSelection, QueryContext } from '../../mol-model/structure';
import { superposeStructures } from '../../mol-model/structure/structure/util/superposition';
import Expression from '../../mol-script/language/expression';

export type SuperpositionTestInput = {
    pdbId: string,
    auth_asym_id: string,
    matrix: Mat4
}[];

// function getAxisAngleTranslation(m: Mat4) {
//     const translation = Mat4.getTranslation(Vec3.zero(), m);
//     const axis = Vec3.zero();
//     const angle = 180 / Math.PI * Quat.getAxisAngle(axis, Mat4.getRotation(Quat.zero(), m));
//     return { translation, axis, angle };
// }

export function buildStaticSuperposition(ctx: PluginContext, src: SuperpositionTestInput) {
    const b = ctx.state.data.build().toRoot();
    for (const s of src) {
        StateHelper.visual(ctx,
            StateHelper.transform(
                StateHelper.selectChain(
                    StateHelper.structure(
                        StateHelper.getModel(StateHelper.download(b, `https://www.ebi.ac.uk/pdbe/static/entry/${s.pdbId}_updated.cif`), 'cif')),
                    s.auth_asym_id
                ),
                s.matrix
            )
        );
    }
    return b;
}

export const StaticSuperpositionTestData: SuperpositionTestInput = [
    { pdbId: '1aj5', auth_asym_id: 'A', matrix: Mat4.identity() },
    { pdbId: '1df0', auth_asym_id: 'B', matrix: Mat4.ofRows([
        [0.406, 0.879, 0.248, -200.633],
        [0.693, -0.473, 0.544, 73.403],
        [0.596, -0.049, -0.802, -14.209],
        [0, 0, 0, 1]] )},
    { pdbId: '1dvi', auth_asym_id: 'A', matrix: Mat4.ofRows([
        [-0.053, -0.077, 0.996, -45.633],
        [-0.312, 0.949, 0.057, -12.255],
        [-0.949, -0.307, -0.074, 53.562],
        [0, 0, 0, 1]] )}
];

export async function dynamicSuperpositionTest(ctx: PluginContext, src: string[], comp_id: string) {
    const state = ctx.state.data;

    const structures = state.build().toRoot();
    for (const s of src) {
        StateHelper.structure(
            StateHelper.getModel(StateHelper.download(structures, `https://www.ebi.ac.uk/pdbe/static/entry/${s}_updated.cif`), 'cif'));
    }

    await PluginCommands.State.Update(ctx, { state, tree: structures });

    const pivot = MS.struct.filter.first([
        MS.struct.generator.atomGroups({
            'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), comp_id]),
            'group-by': MS.struct.atomProperty.macromolecular.residueKey()
        })
    ]);
    const rest = MS.struct.modifier.exceptBy({
        0: MS.struct.generator.all(),
        by: pivot
    });

    const query = compile<StructureSelection>(pivot);
    const xs = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure));
    const selections = xs.map(s => StructureSelection.toLociWithCurrentUnits(query(new QueryContext(s.obj!.data))));

    const transforms = superposeStructures(selections);
    const visuals = state.build();

    siteVisual(ctx, StateHelper.selectSurroundingsOfFirstResidue(visuals.to(xs[0].transform.ref), 'HEM', 7), pivot, rest);
    for (let i = 1; i < selections.length; i++) {
        const root = visuals.to(xs[i].transform.ref);
        siteVisual(ctx,
            StateHelper.transform(StateHelper.selectSurroundingsOfFirstResidue(root, 'HEM', 7), transforms[i - 1].bTransform),
            pivot, rest);
    }

    await PluginCommands.State.Update(ctx, { state, tree: visuals });
}

function siteVisual(ctx: PluginContext, b: StateBuilder.To<PSO.Molecule.Structure>, pivot: Expression, rest: Expression) {
    StateHelper.ballsAndSticks(ctx, b, pivot, 'residue-name');
    StateHelper.ballsAndSticks(ctx, b, rest, 'uniform');
}