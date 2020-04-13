/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder, ComplexRepresentation } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { LabelTextVisual, LabelTextParams } from '../visual/label-text';
import { MarkerAction } from '../../../mol-util/marker-action';

const LabelVisuals = {
    'label-text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, LabelTextParams>) => ComplexRepresentation('Label text', ctx, getParams, LabelTextVisual),
};

export const LabelParams = {
    ...LabelTextParams,
    visuals: PD.MultiSelect(['label-text'], PD.objectToOptions(LabelVisuals)),
};
export type LabelParams = typeof LabelParams
export function getLabelParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(LabelParams);
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
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});