/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { TextBuilder } from '@molstar/graphics/geo/geometry/text/text-builder';
import { Structure } from '@molstar/model/model/structure';
import { ComplexTextVisual, ComplexVisual } from '@molstar/graphics/repr/structure/complex-visual';
import * as Original from '@molstar/graphics/repr/structure/visual/label-text';
import {
  eachSerialElement,
  ElementIterator,
  getSerialElementLoci,
} from '@molstar/graphics/repr/structure/visual/util/element';
import { VisualUpdateState } from '@molstar/graphics/repr/util';
import type { VisualContext } from '@molstar/graphics/repr/visual';
import { Theme } from '@molstar/graphics/theme/theme';
import { arrayEqual } from '@molstar/core/util';
import { ColorNames } from '@molstar/core/util/color/names';
import { omitObjectKeys } from '@molstar/core/util/object';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { FormatTemplate } from '@molstar/core/util/string-format';
import { textPropsForSelection } from '@molstar/mvs/helpers/label-text';
import type { MVSAnnotationRow } from '@molstar/mvs/helpers/schemas';
import { GroupedArray } from '@molstar/mvs/helpers/utils';
import { getMVSAnnotationForStructure, MVSAnnotation } from '../annotation-prop.js';

/** Parameter definition for "label-text" visual in "MVS Annotation Label" representation */
export type MVSAnnotationLabelTextParams = typeof MVSAnnotationLabelTextParams;
export const MVSAnnotationLabelTextParams = {
  annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property', isEssential: true }),
  fieldName: PD.Text('label', {
    description: 'Annotation field (column) from which to take label contents',
    isEssential: true,
  }),
  textFormat: PD.Text('{}', {
    description:
      'Formatting template for the label text. Supports simplified f-string syntax. May reference multiple annotation fields. If value in any field is not defined, label will not be displayed.',
    isEssential: true,
  }),
  groupByFields: PD.ObjectList({ fieldName: PD.Text() }, (obj) => obj.fieldName, {
    defaultValue: [{ fieldName: 'group_id' }],
    description:
      'Set of annotation fields for grouping annotation rows into label instances (i.e. annotation rows with the same values in all group-by fields will yield one label instance). Annotation row with undefined value in any group-by field is considered a separate label instance.',
    isEssential: true,
  }),
  ...omitObjectKeys(Original.LabelTextParams, ['level', 'chainScale', 'residueScale', 'elementScale']),
  borderColor: { ...Original.LabelTextParams.borderColor, defaultValue: ColorNames.black },
};

/** Parameter values for "label-text" visual in "MVS Annotation Label" representation */
export type MVSAnnotationLabelTextProps = PD.Values<MVSAnnotationLabelTextParams>;

/** Create "label-text" visual for "MVS Annotation Label" representation */
export function MVSAnnotationLabelTextVisual(materialId: number): ComplexVisual<MVSAnnotationLabelTextParams> {
  return ComplexTextVisual<MVSAnnotationLabelTextParams>(
    {
      defaultProps: PD.getDefaultValues(MVSAnnotationLabelTextParams),
      createGeometry: createLabelText,
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<MVSAnnotationLabelTextParams>,
        currentProps: PD.Values<MVSAnnotationLabelTextParams>,
      ) => {
        state.createGeometry =
          newProps.annotationId !== currentProps.annotationId ||
          newProps.fieldName !== currentProps.fieldName ||
          newProps.textFormat !== currentProps.textFormat ||
          !arrayEqual(newProps.groupByFields, currentProps.groupByFields);
      },
    },
    materialId,
  );
}

function createLabelText(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: MVSAnnotationLabelTextProps,
  text?: Text,
): Text {
  const { annotation, model } = getMVSAnnotationForStructure(structure, props.annotationId);
  const rows = annotation?.getRows() ?? [];
  const groups = GroupedArray.groupIndices(
    rows,
    rowGroupingFunction(
      annotation!,
      props.groupByFields.map((x) => x.fieldName),
    ),
  );
  const builder = TextBuilder.create(props, groups.count, groups.count / 2, text);
  const template = FormatTemplate(props.textFormat);
  for (let iGroup = 0; iGroup < groups.count; iGroup++) {
    const rowIndicesInGroup = GroupedArray.getGroup(groups, iGroup);
    const labelText = template.format((field) =>
      annotation!.getValueForRow(rowIndicesInGroup[0], field || props.fieldName),
    );
    if (!labelText) continue;
    const rowsInGroup = rowIndicesInGroup.map((i) => rows[i]);
    const p = textPropsForSelection(structure, rowsInGroup, model);
    if (!p) continue;
    builder.add(labelText, p.center[0], p.center[1], p.center[2], p.depth, p.scale, p.group);
  }
  return builder.getText();
}

function rowGroupingFunction(
  annotation: MVSAnnotation,
  groupByFields: string[],
): (row: MVSAnnotationRow, i: number) => string | undefined {
  if (groupByFields.length === 1) {
    const groupByField = groupByFields[0];
    return (row, i) => annotation.getValueForRow(i, groupByField);
  }
  if (groupByFields.length === 0) {
    return () => '';
  }
  return (row, i) => {
    const values = groupByFields.map((field) => annotation.getValueForRow(i, field));
    if (values.includes(undefined)) return undefined;
    return values.join('\t');
  };
}
