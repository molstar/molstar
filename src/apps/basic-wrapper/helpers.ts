/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { PluginContext } from '../../mol-plugin/context';
import { PluginStateObject as PSO } from '../../mol-plugin/state/objects';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StructureRepresentation3DHelpers } from '../../mol-plugin/state/transforms/representation';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { StateBuilder } from '../../mol-state';
import Expression from '../../mol-script/language/expression';
import { BuiltInColorThemeName } from '../../mol-theme/color';
type SupportedFormats = 'cif' | 'pdb'

export namespace StateHelper {
    export function download(b: StateBuilder.To<PSO.Root>, url: string, ref?: string) {
        return b.apply(StateTransforms.Data.Download, { url, isBinary: false }, { ref });
    }

    export function getModel(b: StateBuilder.To<PSO.Data.Binary | PSO.Data.String>, format: SupportedFormats, modelIndex = 0) {
        const parsed = format === 'cif'
            ? b.apply(StateTransforms.Data.ParseCif).apply(StateTransforms.Model.TrajectoryFromMmCif)
            : b.apply(StateTransforms.Model.TrajectoryFromPDB);

        return parsed.apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex });
    }

    export function structure(b: StateBuilder.To<PSO.Molecule.Model>) {
        return b.apply(StateTransforms.Model.StructureFromModel, void 0, { tags: 'structure' })
    };

    export function selectChain(b: StateBuilder.To<PSO.Molecule.Structure>, auth_asym_id: string) {
        const expression = MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), auth_asym_id])
        })
        return b.apply(StateTransforms.Model.StructureSelectionFromExpression, { expression, label: `Chain ${auth_asym_id}` });
    }

    export function select(b: StateBuilder.To<PSO.Molecule.Structure>, expression: Expression) {
        return b.apply(StateTransforms.Model.StructureSelectionFromExpression, { expression });
    }

    export function selectSurroundingsOfFirstResidue(b: StateBuilder.To<PSO.Molecule.Structure>, comp_id: string, radius: number) {
        const expression = MS.struct.modifier.includeSurroundings({
            0: MS.struct.filter.first([
                MS.struct.generator.atomGroups({
                    'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), comp_id]),
                    'group-by': MS.struct.atomProperty.macromolecular.residueKey()
                })
            ]),
            radius
        })
        return b.apply(StateTransforms.Model.StructureSelectionFromExpression, { expression, label: `Surr. ${comp_id} (${radius} ang)` });
    }

    export function identityTransform(b: StateBuilder.To<PSO.Molecule.Structure>, m: Mat4) {
        return b.apply(StateTransforms.Model.TransformStructureConformation,
            { axis: Vec3.create(1, 0, 0), angle: 0, translation: Vec3.zero() },
            { tags: 'transform' });
    }

    export function transform(b: StateBuilder.To<PSO.Molecule.Structure>, matrix: Mat4) {
        return b.apply(StateTransforms.Model.TransformStructureConformationByMatrix, { matrix }, { tags: 'transform' });
    }

    export function assemble(b: StateBuilder.To<PSO.Molecule.Model>, id?: string) {
        return b.apply(StateTransforms.Model.StructureAssemblyFromModel, { id: id || 'deposited' }, { tags: 'asm' })
    }

    export function visual(ctx: PluginContext, visualRoot: StateBuilder.To<PSO.Molecule.Structure>) {
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'cartoon'), { tags: 'seq-visual' });
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick'), { tags: 'het-visual' });
        // visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
        //     .apply(StateTransforms.Representation.StructureRepresentation3D,
        //         StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick', { alpha: 0.51 }), { tags: 'water-visual' });
        return visualRoot;
    }

    export function ballsAndSticks(ctx: PluginContext, visualRoot: StateBuilder.To<PSO.Molecule.Structure>, expression: Expression, coloring?: BuiltInColorThemeName) {
        visualRoot
            .apply(StateTransforms.Model.StructureSelectionFromExpression, { expression })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick', void 0, coloring), { tags: 'het-visual' });
        return visualRoot;
    }

}