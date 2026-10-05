/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Representation, type RepresentationContext, type RepresentationParamsGetter } from '@molstar/graphics/repr/representation';
import { Structure } from '@molstar/model/model/structure';
import { type StructureRepresentation, StructureRepresentationStateBuilder } from '@molstar/graphics/repr/structure/representation';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';

// TODO: Convert this approach to a more usual one that creates visuals during parse and shows them
// during visuals.

const KinemageDataVisuals = {
};

export const KinemageDataParams = {
    visuals: PD.MultiSelect([], PD.objectToOptions(KinemageDataVisuals)),
};
export type KinemageDataParams = typeof KinemageDataParams
export type KinemageDataProps = PD.Values<KinemageDataParams>

export function getKinemageDataParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(KinemageDataParams);
}

export type KinemageDataRepresentation = StructureRepresentation<KinemageDataParams>
export function KinemageDataRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, KinemageDataParams>): KinemageDataRepresentation {
    return Representation.createMulti('Membrane Orientation', ctx, getParams, StructureRepresentationStateBuilder, KinemageDataVisuals as unknown as Representation.Def<Structure, KinemageDataParams>);
}
