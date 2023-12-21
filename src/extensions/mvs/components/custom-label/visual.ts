/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { SortedArray } from '../../../../mol-data/int';
import { Text } from '../../../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../../../mol-geo/geometry/text/text-builder';
import { Structure } from '../../../../mol-model/structure';
import { ComplexTextVisual, ComplexVisual } from '../../../../mol-repr/structure/complex-visual';
import * as Original from '../../../../mol-repr/structure/visual/label-text';
import { ElementIterator, eachSerialElement, getSerialElementLoci } from '../../../../mol-repr/structure/visual/util/element';
import { VisualUpdateState } from '../../../../mol-repr/util';
import { VisualContext } from '../../../../mol-repr/visual';
import { Theme } from '../../../../mol-theme/theme';
import { deepEqual } from '../../../../mol-util';
import { ColorNames } from '../../../../mol-util/color/names';
import { omitObjectKeys } from '../../../../mol-util/object';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { textPropsForSelection } from '../../helpers/label-text';
import { SelectorParams, substructureFromSelector } from '../selector';


/** Parameter definition for "label-text" visual in "Custom Label" representation */
export type CustomLabelTextParams = typeof CustomLabelTextParams
export const CustomLabelTextParams = {
    items: PD.ObjectList(
        {
            text: PD.Text('¯\\_(ツ)_/¯'),
            position: PD.MappedStatic('selection', {
                x_y_z: PD.Group({
                    x: PD.Numeric(0),
                    y: PD.Numeric(0),
                    z: PD.Numeric(0),
                    scale: PD.Numeric(1, { min: 0, max: 20, step: 0.1 })
                }),
                selection: PD.Group({
                    selector: SelectorParams,
                }),
            }),
        },
        obj => obj.text,
        { isEssential: true }
    ),
    ...omitObjectKeys(Original.LabelTextParams, ['level', 'chainScale', 'residueScale', 'elementScale']),
    borderColor: { ...Original.LabelTextParams.borderColor, defaultValue: ColorNames.black },
};

/** Parameter values for "label-text" visual in "Custom Label" representation */
export type CustomLabelTextProps = PD.Values<CustomLabelTextParams>

/** Create "label-text" visual for "Custom Label" representation */
export function CustomLabelTextVisual(materialId: number): ComplexVisual<CustomLabelTextParams> {
    return ComplexTextVisual<CustomLabelTextParams>({
        defaultProps: PD.getDefaultValues(CustomLabelTextParams),
        createGeometry: createLabelText,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CustomLabelTextParams>, currentProps: PD.Values<CustomLabelTextParams>) => {
            state.createGeometry = !deepEqual(newProps.items, currentProps.items);
        }
    }, materialId);
}

function createLabelText(ctx: VisualContext, structure: Structure, theme: Theme, props: CustomLabelTextProps, text?: Text): Text {
    const count = props.items.length;
    const builder = TextBuilder.create(props, count, count / 2, text);
    for (const item of props.items) {
        switch (item.position.name) {
            case 'x_y_z':
                const scale = item.position.params.scale;
                builder.add(item.text, item.position.params.x, item.position.params.y, item.position.params.z, scale, scale, 0);
                break;
            case 'selection':
                const substructure = substructureFromSelector(structure, item.position.params.selector);
                const p = textPropsForSelection(substructure, theme.size.size, {});
                const group = serialIndexOfSubstructure(structure, substructure) ?? 0;
                if (p) builder.add(item.text, p.center[0], p.center[1], p.center[2], p.depth, p.scale, group);
                break;
        }
    }
    return builder.getText();
}

/** Return the serial index within `structure` of the first element of `substructure` (or `undefined` in that element is not in `structure`)  */
function serialIndexOfSubstructure(structure: Structure, substructure: Structure): number | undefined {
    if (substructure.isEmpty) return undefined;
    const theUnit = substructure.units[0];
    const theElement = theUnit.elements[0];
    for (const unit of structure.units) {
        if (unit.model.id === theUnit.model.id && SortedArray.has(unit.elements, theElement)) {
            return structure.serialMapping.getSerialIndex(unit, theElement);
        }
    }
    return undefined;
}
