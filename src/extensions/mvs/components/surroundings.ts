/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { StructureElement } from '../../../mol-model/structure';
import { createStructureComponent } from '../../../mol-plugin-state/helpers/structure-component';
import { PluginStateTransform, PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';


export const StructureSurroundingsParams = {
    radius: PD.Numeric(5, { min: 0 }, { description: 'Surroundings radius in Angstroms' }),
    includeSelf: PD.Boolean(true, { description: 'Include parent selection itself in the surroundings' }),
    wholeResidues: PD.Boolean(true, { description: 'Include whole residues, instead of individual atoms' }),
    nullIfEmpty: PD.Optional(PD.Boolean(true, { isHidden: true })),
};
export type StructureSurroundingsParams = typeof StructureSurroundingsParams;
export type StructureSurroundingsProps = PD.ValuesFor<StructureSurroundingsParams>;


export type StructureSurroundings = typeof StructureSurroundings;
export const StructureSurroundings = PluginStateTransform.BuiltIn({
    name: 'structure-surroundings',
    display: { name: 'Surroundings', description: 'Surroundings of a structure component.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: StructureSurroundingsParams,
})({
    apply({ a, params, cache }) {
        const struct = a.data;
        const rootStruct = struct.parent ?? struct;
        const targetBundle = StructureElement.Bundle.fromSubStructure(rootStruct, struct);
        const targetExpr = StructureElement.Bundle.toExpression(targetBundle);
        let surroundingsExpr = MolScriptBuilder.struct.modifier.includeSurroundings({
            0: targetExpr,
            radius: params.radius,
            'as-whole-residues': params.wholeResidues,
        });
        if (!params.includeSelf) {
            surroundingsExpr = MolScriptBuilder.struct.modifier.exceptBy({
                0: surroundingsExpr,
                by: targetExpr,
            });
        }
        return createStructureComponent(rootStruct, { label: `Surroundings (${params.radius} Å)`, type: { name: 'expression', params: surroundingsExpr }, nullIfEmpty: params.nullIfEmpty }, cache as any);
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});
