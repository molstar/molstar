/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../context';
import { StateBuilder } from '../../mol-state';
import { PluginStateObject } from '../state/objects';
import { StateTransforms } from '../state/transforms';
import { StructureRepresentation3DHelpers } from '../state/transforms/representation';
import { StructureComplexElementTypes } from '../state/transforms/model';

export function createDefaultStructureComplex(
    ctx: PluginContext, root: StateBuilder.To<PluginStateObject.Molecule.Structure>
) {
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'protein-or-nucleic' }, { tags: StructureComplexElementTypes['protein-or-nucleic'] })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'cartoon'));

    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'ligand' }, { tags: StructureComplexElementTypes.ligand })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick'));

    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'modified' }, { tags: StructureComplexElementTypes.modified })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick', void 0, 'polymer-id'));

    const branched = root.apply(StateTransforms.Model.StructureComplexElement, { type: 'branched' }, { tags: StructureComplexElementTypes.branched })

    branched.apply(StateTransforms.Representation.StructureRepresentation3D,
        StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick', { alpha: 0.15 }));
    branched.apply(StateTransforms.Representation.StructureRepresentation3D,
        StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'carbohydrate'));

    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' }, { tags: StructureComplexElementTypes.water })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick', { alpha: 0.51 }));

    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'coarse' }, { tags: StructureComplexElementTypes.coarse })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'spacefill', {}, 'polymer-id'));
}