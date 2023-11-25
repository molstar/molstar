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
import { deepEqual } from '../../../../mol-util';
import { ColorNames } from '../../../../mol-util/color/names';
import { omitObjectKeys } from '../../../../mol-util/object';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { textPropsForSelection } from '../../helpers/label-text';
import { MaybeIntegerParamDefinition, MaybeStringParamDefinition } from '../../helpers/param-definition';


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
                    label_entity_id: MaybeStringParamDefinition(),
                    label_asym_id: MaybeStringParamDefinition(),
                    auth_asym_id: MaybeStringParamDefinition(),

                    label_seq_id: MaybeIntegerParamDefinition(),
                    auth_seq_id: MaybeIntegerParamDefinition(),
                    pdbx_PDB_ins_code: MaybeStringParamDefinition(),
                    /** Minimum label_seq_id (inclusive) */
                    beg_label_seq_id: MaybeIntegerParamDefinition(undefined, { description: 'Minimum label_seq_id (inclusive)' }),
                    /** Maximum label_seq_id (inclusive) */
                    end_label_seq_id: MaybeIntegerParamDefinition(),
                    /** Minimum auth_seq_id (inclusive) */
                    beg_auth_seq_id: MaybeIntegerParamDefinition(),
                    /** Maximum auth_seq_id (inclusive) */
                    end_auth_seq_id: MaybeIntegerParamDefinition(),

                    /** Atom name like 'CA', 'N', 'O'... */
                    label_atom_id: MaybeStringParamDefinition(),
                    /** Atom name like 'CA', 'N', 'O'... */
                    auth_atom_id: MaybeStringParamDefinition(),
                    /** Element symbol like 'H', 'HE', 'LI', 'BE'... */
                    type_symbol: MaybeStringParamDefinition(),
                    /** Unique atom identifier across conformations (_atom_site.id) */
                    atom_id: MaybeIntegerParamDefinition(),
                    /** 0-base index of the atom in the source data */
                    atom_index: MaybeIntegerParamDefinition(),

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
                const p = textPropsForSelection(structure, theme.size.size, item.position.params);
                if (p) builder.add(item.text, p.center[0], p.center[1], p.center[2], p.depth, p.scale, p.group);
                break;
        }
    }
    return builder.getText();
}
