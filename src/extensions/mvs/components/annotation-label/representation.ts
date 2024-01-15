/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Structure } from '../../../../mol-model/structure';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../../mol-repr/representation';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../../../../mol-repr/structure/representation';
import { MarkerAction } from '../../../../mol-util/marker-action';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { isMVSStructure } from '../is-mvs-model-prop';
import { MVSAnnotationLabelTextParams, MVSAnnotationLabelTextVisual } from './visual';


/** Components of "MVS Annotation Label" representation */
const MVSAnnotationLabelVisuals = {
    'label-text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MVSAnnotationLabelTextParams>) => ComplexRepresentation('Label text', ctx, getParams, MVSAnnotationLabelTextVisual),
};

/** Parameter definition for representation type "MVS Annotation Label" */
export type MVSAnnotationLabelParams = typeof MVSAnnotationLabelParams
export const MVSAnnotationLabelParams = {
    ...MVSAnnotationLabelTextParams,
    visuals: PD.MultiSelect(['label-text'], PD.objectToOptions(MVSAnnotationLabelVisuals)),
};

/** Parameter values for representation type "MVS Annotation Label" */
export type MVSAnnotationLabelProps = PD.ValuesFor<MVSAnnotationLabelParams>

/** Structure representation type "MVS Annotation Label", allowing showing labels based on "MVS Annotations" custom props */
export type MVSAnnotationLabelRepresentation = StructureRepresentation<MVSAnnotationLabelParams>
export function MVSAnnotationLabelRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MVSAnnotationLabelParams>): MVSAnnotationLabelRepresentation {
    const repr = Representation.createMulti('Label', ctx, getParams, StructureRepresentationStateBuilder, MVSAnnotationLabelVisuals as unknown as Representation.Def<Structure, MVSAnnotationLabelParams>);
    repr.setState({ pickable: false, markerActions: MarkerAction.None });
    return repr;
}

/** A thingy that is needed to register representation type "MVS Annotation Label", allowing showing labels based on "MVS Annotations" custom props */
export const MVSAnnotationLabelRepresentationProvider = StructureRepresentationProvider({
    name: 'mvs-annotation-label',
    label: 'MVS Annotation Label',
    description: 'Displays labels based on annotation custom model property',
    factory: MVSAnnotationLabelRepresentation,
    getParams: () => MVSAnnotationLabelParams,
    defaultValues: PD.getDefaultValues(MVSAnnotationLabelParams),
    defaultColorTheme: { name: 'uniform' }, // this ain't workin
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0 && isMVSStructure(structure),
});
