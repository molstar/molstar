/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { InteractionType } from '../../mol-model-props/computed/interactions/common';
import { StructureElement } from '../../mol-model/structure';
import { StructureElementSchema, StructureElementSchemaItem } from '../../mol-model/structure/query/schema';

interface InteractionElementSchemaBase {
    aStructureRef?: string,
    a: StructureElementSchema,
    bStructureRef?: string,
    b: StructureElementSchema,
    description?: string,
}

export type InteractionElementSchema =
    | { kind: 'unknown' } & InteractionElementSchemaBase
    | { kind: 'ionic' } & InteractionElementSchemaBase
    | { kind: 'pi-stacking' } & InteractionElementSchemaBase
    | { kind: 'cation-pi' } & InteractionElementSchemaBase
    | { kind: 'halogen-bond' } & InteractionElementSchemaBase
    | { kind: 'hydrogen-bond', hydrogenStructureRef?: string, hydrogen?: StructureElementSchemaItem } & InteractionElementSchemaBase
    | { kind: 'weak-hydrogen-bond', hydrogenStructureRef?: string, hydrogen?: StructureElementSchemaItem } & InteractionElementSchemaBase
    | { kind: 'hydrophobic' } & InteractionElementSchemaBase
    | { kind: 'metal-coordination' } & InteractionElementSchemaBase
    | { kind: 'salt-bridge' } & InteractionElementSchemaBase
    | { kind: 'covalent', degree?: number, aromatic?: boolean } & InteractionElementSchemaBase

export type InteractionKind = InteractionElementSchema['kind']

export const InteractionKinds: InteractionKind[] = [
    'unknown',
    'ionic',
    'pi-stacking',
    'cation-pi',
    'halogen-bond',
    'hydrogen-bond',
    'weak-hydrogen-bond',
    'hydrophobic',
    'metal-coordination',
    'salt-bridge',
    'covalent',
];

export type InteractionInfo =
    | { kind: 'unknown' }
    | { kind: 'ionic' }
    | { kind: 'pi-stacking' }
    | { kind: 'cation-pi' }
    | { kind: 'halogen-bond' }
    | { kind: 'hydrogen-bond', hydrogenStructureRef?: string, hydrogen?: StructureElement.Loci }
    | { kind: 'weak-hydrogen-bond', hydrogenStructureRef?: string, hydrogen?: StructureElement.Loci }
    | { kind: 'hydrophobic' }
    | { kind: 'metal-coordination' }
    | { kind: 'salt-bridge' }
    | { kind: 'covalent', degree?: number, aromatic?: boolean }

export interface StructureInteractionElement {
    // Pass the schema when loading from custom data
    sourceSchema?: InteractionElementSchema,

    info: InteractionInfo,
    aStructureRef?: string,
    a: StructureElement.Loci,
    bStructureRef?: string,
    b: StructureElement.Loci,
}

export interface StructureInteractions {
    kind: 'structure-interactions',
    elements: StructureInteractionElement[],
}

export const InteractionTypeToKind = {
    [InteractionType.Unknown]: 'unknown' as InteractionKind,
    [InteractionType.Ionic]: 'ionic' as InteractionKind,
    [InteractionType.CationPi]: 'cation-pi' as InteractionKind,
    [InteractionType.PiStacking]: 'pi-stacking' as InteractionKind,
    [InteractionType.HydrogenBond]: 'hydrogen-bond' as InteractionKind,
    [InteractionType.HalogenBond]: 'halogen-bond' as InteractionKind,
    [InteractionType.Hydrophobic]: 'hydrophobic' as InteractionKind,
    [InteractionType.MetalCoordination]: 'metal-coordination' as InteractionKind,
    [InteractionType.WeakHydrogenBond]: 'weak-hydrogen-bond' as InteractionKind,
};