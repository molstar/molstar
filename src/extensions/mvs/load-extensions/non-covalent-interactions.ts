/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { StructureSurroundings } from '../components/surroundings';
import { MolstarLoadingExtension } from '../load';
import { UpdateTarget } from '../load-generic';
import { getCustomProps } from '../tree/generic/tree-schema';


const DefaultNonCovalentInteractionRadius = 5;

export const NonCovalentInteractionsExtension: MolstarLoadingExtension<{}> = {
    id: 'wwpdb/non-covalent-interactions',
    description: 'Allow showing non-covalent interactions around components with molstar_show_non_covalent_interactions additional property',
    createExtensionContext: () => ({}),
    action: (updateTarget, node, context, extContext) => {
        if (node.kind !== 'component' && node.kind !== 'component_from_uri' && node.kind !== 'component_from_source') return;

        type CustomProps = {
            molstar_show_non_covalent_interactions?: boolean,
            molstar_non_covalent_interactions_radius_ang?: number,
        };
        const customProps = getCustomProps<CustomProps>(node);
        if (!customProps.molstar_show_non_covalent_interactions) return undefined;

        const surroundings = UpdateTarget.apply(updateTarget, StructureSurroundings, {
            radius: customProps.molstar_non_covalent_interactions_radius_ang ?? DefaultNonCovalentInteractionRadius,
            includeSelf: true,
            wholeResidues: true,
            nullIfEmpty: false,
        });
        // Bubble on target
        UpdateTarget.apply(updateTarget, StateTransforms.Representation.StructureRepresentation3D, {
            type: { name: 'ball-and-stick', params: { sizeFactor: 0.22, sizeAspectRatio: 0.73, adjustCylinderLength: true, xrayShaded: true, aromaticBonds: false, multipleBonds: 'off', excludeTypes: ['hydrogen-bond', 'metal-coordination'] } },
            colorTheme: { name: 'element-symbol', params: {} },
            sizeTheme: { name: 'physical', params: {} },
        });
        // Ball-and-stick on surrounding
        UpdateTarget.apply(surroundings, StateTransforms.Representation.StructureRepresentation3D, {
            type: { name: 'ball-and-stick', params: { sizeFactor: 0.16, excludeTypes: ['hydrogen-bond', 'metal-coordination'] } },
            colorTheme: { name: 'element-symbol', params: {} },
            sizeTheme: { name: 'physical', params: {} },
        });
        // Non-covalent interactions
        UpdateTarget.apply(surroundings, StateTransforms.Representation.StructureRepresentation3D, {
            type: { name: 'interactions', params: {} },
            colorTheme: { name: 'interaction-type', params: {} },
            sizeTheme: { name: 'uniform', params: {} },
        });
    },
};
