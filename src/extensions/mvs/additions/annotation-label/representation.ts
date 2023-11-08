/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Structure } from '../../../../mol-model/structure';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../../mol-repr/representation';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../../../../mol-repr/structure/representation';
import { MarkerAction } from '../../../../mol-util/marker-action';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { AnnotationLabelTextParams, AnnotationLabelTextVisual } from './visual';


/** Components of "Annotation Label" representation */
const AnnotationLabelVisuals = {
    'label-text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, AnnotationLabelTextParams>) => ComplexRepresentation('Label text', ctx, getParams, AnnotationLabelTextVisual),
};

/** Parameter definition for representation type "Annotation Label" */
export type AnnotationLabelParams = typeof AnnotationLabelParams
export const AnnotationLabelParams = {
    ...AnnotationLabelTextParams,
    visuals: PD.MultiSelect(['label-text'], PD.objectToOptions(AnnotationLabelVisuals)),
};

/** Parameter values for representation type "Annotation Label" */
export type AnnotationLabelProps = PD.ValuesFor<AnnotationLabelParams>

/** Structure representation type "Annotation Label", allowing showing labels based on "Annotations" custom props */
export type AnnotationLabelRepresentation = StructureRepresentation<AnnotationLabelParams>
export function AnnotationLabelRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, AnnotationLabelParams>): AnnotationLabelRepresentation {
    const repr = Representation.createMulti('Label', ctx, getParams, StructureRepresentationStateBuilder, AnnotationLabelVisuals as unknown as Representation.Def<Structure, AnnotationLabelParams>);
    repr.setState({ pickable: false, markerActions: MarkerAction.None });
    return repr;
}

/** A thingy that is needed to register representation type "Annotation Label", allowing showing labels based on "Annotations" custom props */
export const AnnotationLabelRepresentationProvider = StructureRepresentationProvider({
    name: 'mvs-annotation-label',
    label: 'Annotation Label',
    description: 'Displays labels based on annotation custom model property',
    factory: AnnotationLabelRepresentation,
    getParams: () => AnnotationLabelParams,
    defaultValues: PD.getDefaultValues(AnnotationLabelParams),
    defaultColorTheme: { name: 'uniform' }, // this ain't workin
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
});
