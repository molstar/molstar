/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { InteractionType } from '../../mol-model-props/computed/interactions/common';
import { StructureElementSchema, StructureElementSchemaItem } from '../../mol-model/structure/query/schema';

interface InteractionElementsSchema {
    aStructureRef?: string,
    a: StructureElementSchema,
    bStructureRef?: string,
    b: StructureElementSchema,
}

export type InteractionSchema =
    | { kind: 'unknown' } & InteractionElementsSchema
    | { kind: 'ionic' } & InteractionElementsSchema
    | { kind: 'pi-stacking' } & InteractionElementsSchema
    | { kind: 'cation-pi' } & InteractionElementsSchema
    | { kind: 'halogen-bond' } & InteractionElementsSchema
    | { kind: 'hydrogen-bond', hydrogen?: StructureElementSchemaItem } & InteractionElementsSchema
    | { kind: 'weak-hydrogen-bond', hydrogen?: StructureElementSchemaItem } & InteractionElementsSchema
    | { kind: 'hydrophobic' } & InteractionElementsSchema
    | { kind: 'metal-coordination' } & InteractionElementsSchema
    | { kind: 'covalent', degree?: number } & InteractionElementsSchema

export type InteractionKind = InteractionSchema['kind']

export interface StructureInteractions {
    elements: InteractionSchema[],
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