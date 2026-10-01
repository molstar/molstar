/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder, ComplexRepresentation } from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { LabelTextVisual, LabelTextParams } from '../visual/label-text.js';
import { MarkerAction } from '@molstar/core/util/marker-action';

const LabelVisuals = {
    'label-text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, LabelTextParams>) => ComplexRepresentation('Label text', ctx, getParams, LabelTextVisual),
};

export const LabelParams = {
    ...LabelTextParams,
    visuals: PD.MultiSelect(['label-text'], PD.objectToOptions(LabelVisuals)),
};
export type LabelParams = typeof LabelParams
export function getLabelParams(ctx: ThemeRegistryContext, structure: Structure) {
    return LabelParams;
}

export type LabelRepresentation = StructureRepresentation<LabelParams>
export function LabelRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, LabelParams>): LabelRepresentation {
    const repr = Representation.createMulti('Label', ctx, getParams, StructureRepresentationStateBuilder, LabelVisuals as unknown as Representation.Def<Structure, LabelParams>);
    repr.setState({ pickable: false, markerActions: MarkerAction.None });
    return repr;
}

export const LabelRepresentationProvider = StructureRepresentationProvider({
    name: 'label',
    label: 'Label',
    description: 'Displays labels.',
    factory: LabelRepresentation,
    getParams: getLabelParams,
    defaultValues: PD.getDefaultValues(LabelParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});