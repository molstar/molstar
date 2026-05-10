/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../mol-repr/representation';
import { Structure } from '../../mol-model/structure';
import { StructureRepresentation, StructureRepresentationStateBuilder } from '../../mol-repr/structure/representation';
import { ThemeRegistryContext } from '../../mol-theme/theme';

/// @todo Convert this approach to a more usual one that creates visuals during parse and shows them
/// during visuals.

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
