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

const LabelVisuals = {
    'label-text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, LabelTextParams>) => ComplexRepresentation('Label text', ctx, getParams, LabelTextVisual),
}
type LabelVisualName = keyof typeof LabelVisuals
const LabelVisualOptions = Object.keys(LabelVisuals).map(name => [name, name] as [LabelVisualName, string])

export const LabelParams = {
    ...LabelTextParams,
    visuals: PD.MultiSelect<LabelVisualName>(['label-text'], LabelVisualOptions),
}
export type LabelParams = typeof LabelParams
export function getLabelParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(LabelParams)
}

export type LabelRepresentation = StructureRepresentation<LabelParams>
export function LabelRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, LabelParams>): LabelRepresentation {
    return Representation.createMulti('Label', ctx, getParams, StructureRepresentationStateBuilder, LabelVisuals as unknown as Representation.Def<Structure, LabelParams>)
}

export const LabelRepresentationProvider: StructureRepresentationProvider<LabelParams> = {
    label: 'Label',
    description: 'Displays labels.',
    factory: LabelRepresentation,
    getParams: getLabelParams,
    defaultValues: PD.getDefaultValues(LabelParams),
    defaultColorTheme: 'uniform',
    defaultSizeTheme: 'uniform',
    isApplicable: (structure: Structure) => structure.elementCount > 0
}