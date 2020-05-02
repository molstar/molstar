/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4 } from '../../mol-math/linear-algebra';
import { QueryContext, StructureSelection } from '../../mol-model/structure';
import { superpose } from '../../mol-model/structure/structure/util/superposition';
import { PluginStateObject as PSO } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import Expression from '../../mol-script/language/expression';
import { compile } from '../../mol-script/runtime/query/compiler';
import { StateObjectRef } from '../../mol-state';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Asset } from '../../mol-util/assets';

export type SuperpositionTestInput = {
    pdbId: string,
    auth_asym_id: string,
    matrix: Mat4
}[];

export function buildStaticSuperposition(plugin: PluginContext, src: SuperpositionTestInput) {
    return plugin.dataTransaction(async () => {
        for (const s of src) {
            const { structure } = await loadStructure(plugin, `https://www.ebi.ac.uk/pdbe/static/entry/${s.pdbId}_updated.cif`, 'mmcif');
            await transform(plugin, structure, s.matrix);
            const chain = await plugin.builders.structure.tryCreateComponentFromExpression(structure, chainSelection(s.auth_asym_id), `Chain ${s.auth_asym_id}`);
            if (chain) await plugin.builders.structure.representation.addRepresentation(chain, { type: 'cartoon' });
        }
    });
}

export const StaticSuperpositionTestData: SuperpositionTestInput = [
    {
        pdbId: '1aj5', auth_asym_id: 'A', matrix: Mat4.identity()
    },
    {
        pdbId: '1df0', auth_asym_id: 'B', matrix: Mat4.ofRows([
            [0.406, 0.879, 0.248, -200.633],
            [0.693, -0.473, 0.544, 73.403],
            [0.596, -0.049, -0.802, -14.209],
            [0, 0, 0, 1]])
    },
    {
        pdbId: '1dvi', auth_asym_id: 'A', matrix: Mat4.ofRows([
            [-0.053, -0.077, 0.996, -45.633],
            [-0.312, 0.949, 0.057, -12.255],
            [-0.949, -0.307, -0.074, 53.562],
            [0, 0, 0, 1]])
    }
];

export function dynamicSuperpositionTest(plugin: PluginContext, src: string[], comp_id: string) {
    return plugin.dataTransaction(async () => {
        for (const s of src) {
            await loadStructure(plugin, `https://www.ebi.ac.uk/pdbe/static/entry/${s}_updated.cif`, 'mmcif');
        }

        const pivot = MS.struct.filter.first([
            MS.struct.generator.atomGroups({
                'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), comp_id]),
                'group-by': MS.struct.atomProperty.macromolecular.residueKey()
            })
        ]);

        const rest = MS.struct.modifier.exceptBy({
            0: MS.struct.modifier.includeSurroundings({
                0: pivot,
                radius: 5
            }),
            by: pivot
        });

        const query = compile<StructureSelection>(pivot);
        const xs = plugin.managers.structure.hierarchy.current.structures;
        const selections = xs.map(s => StructureSelection.toLociWithCurrentUnits(query(new QueryContext(s.cell.obj!.data))));

        const transforms = superpose(selections);

        await siteVisual(plugin, xs[0].cell, pivot, rest);
        for (let i = 1; i < selections.length; i++) {
            await transform(plugin, xs[i].cell, transforms[i - 1].bTransform);
            await siteVisual(plugin, xs[i].cell, pivot, rest);
        }
    });
}

async function siteVisual(plugin: PluginContext, s: StateObjectRef<PSO.Molecule.Structure>, pivot: Expression, rest: Expression) {
    const center = await plugin.builders.structure.tryCreateComponentFromExpression(s, pivot, 'pivot');
    if (center) await plugin.builders.structure.representation.addRepresentation(center, { type: 'ball-and-stick', color: 'residue-name' });

    const surr = await plugin.builders.structure.tryCreateComponentFromExpression(s, rest, 'rest');
    if (surr) await plugin.builders.structure.representation.addRepresentation(surr, { type: 'ball-and-stick', color: 'uniform', size: 'uniform', sizeParams: { value: 0.33 } });
}

async function loadStructure(plugin: PluginContext, url: string, format: BuiltInTrajectoryFormat, assemblyId?: string) {
    const data = await plugin.builders.data.download({ url: Asset.Url(url) });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model, assemblyId ? { name: 'assembly', params: { id: assemblyId } } : void 0);

    return { data, trajectory, model, structure };
}

function chainSelection(auth_asym_id: string) {
    return MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), auth_asym_id])
    });
}

function transform(plugin: PluginContext, s: StateObjectRef<PSO.Molecule.Structure>, matrix: Mat4) {
    const b = plugin.state.data.build().to(s)
        .insert(StateTransforms.Model.TransformStructureConformation, { transform: { name: 'matrix', params: { data: matrix, transpose: false } } });
    return plugin.runTask(plugin.state.data.updateTree(b));
}