/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Text } from '../../../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../../../mol-geo/geometry/text/text-builder';
import { Structure } from '../../../../mol-model/structure';
import { ComplexTextVisual, ComplexVisual } from '../../../../mol-repr/structure/complex-visual';
import * as Original from '../../../../mol-repr/structure/visual/label-text';
import { ElementIterator, eachSerialElement, getSerialElementLoci } from '../../../../mol-repr/structure/visual/util/element';
import { VisualUpdateState } from '../../../../mol-repr/util';
import { VisualContext } from '../../../../mol-repr/visual';
import { Theme } from '../../../../mol-theme/theme';
import { ColorNames } from '../../../../mol-util/color/names';
import { omitObjectKeys } from '../../../../mol-util/object';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { textPropsForSelection } from '../../helpers/label-text';
import { groupRows } from '../../helpers/selections';
import { getMVSAnnotationForStructure } from '../annotation-prop';


/** Parameter definition for "label-text" visual in "MVS Annotation Label" representation */
export type MVSAnnotationLabelTextParams = typeof MVSAnnotationLabelTextParams
export const MVSAnnotationLabelTextParams = {
    annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property', isEssential: true }),
    fieldName: PD.Text('label', { description: 'Annotation field (column) from which to take label contents', isEssential: true }),
    ...omitObjectKeys(Original.LabelTextParams, ['level', 'chainScale', 'residueScale', 'elementScale']),
    borderColor: { ...Original.LabelTextParams.borderColor, defaultValue: ColorNames.black },
};

/** Parameter values for "label-text" visual in "MVS Annotation Label" representation */
export type MVSAnnotationLabelTextProps = PD.Values<MVSAnnotationLabelTextParams>

/** Create "label-text" visual for "MVS Annotation Label" representation */
export function MVSAnnotationLabelTextVisual(materialId: number): ComplexVisual<MVSAnnotationLabelTextParams> {
    return ComplexTextVisual<MVSAnnotationLabelTextParams>({
        defaultProps: PD.getDefaultValues(MVSAnnotationLabelTextParams),
        createGeometry: createLabelText,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<MVSAnnotationLabelTextParams>, currentProps: PD.Values<MVSAnnotationLabelTextParams>) => {
            state.createGeometry = newProps.annotationId !== currentProps.annotationId || newProps.fieldName !== currentProps.fieldName;
        }
    }, materialId);
}

function createLabelText(ctx: VisualContext, structure: Structure, theme: Theme, props: MVSAnnotationLabelTextProps, text?: Text): Text {
    const { annotation, model } = getMVSAnnotationForStructure(structure, props.annotationId);
    const rows = annotation?.getRows() ?? [];
    const { count, offsets, grouped } = groupRows(rows);
    const builder = TextBuilder.create(props, count, count / 2, text);
    for (let iGroup = 0; iGroup < count; iGroup++) {
        const iFirstRowInGroup = grouped[offsets[iGroup]];
        const labelText = annotation!.getValueForRow(iFirstRowInGroup, props.fieldName);
        if (!labelText) continue;
        const rowsInGroup = grouped.slice(offsets[iGroup], offsets[iGroup + 1]).map(j => rows[j]);
        const p = textPropsForSelection(structure, theme.size.size, rowsInGroup, model);
        if (!p) continue;
        builder.add(labelText, p.center[0], p.center[1], p.center[2], p.depth, p.scale, p.group);
    }
    return builder.getText();
}
