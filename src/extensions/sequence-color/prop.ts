/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { ElementIndex, Structure } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ElementSet, Selector, SelectorParams } from '../mvs/components/selector';


export namespace SequenceColorProperty {
    /** Parameter definition for custom structure property `SequenceColorProperty` */
    export type Params = typeof Params
    export const Params = {
        colors: PD.ObjectList(
            {
                color: PD.Color(ColorNames.grey, { description: 'Color to apply to a substructure' }),
                selector: SelectorParams,
            },
            obj => Color.toHexStyle(obj.color),
            { description: 'List of substructure-color assignments' }
        ),
    };

    /** Parameter values of custom structure property `SequenceColorProperty` */
    export type Props = PD.Values<Params>

    /** Values of custom structure property `SequenceColorProperty` */
    export interface Data {
        items: {
            selector: Selector,
            color: Color,
            elementSet?: ElementSet,
        }[],
        colorCache: {
            [unitId: number]: {
                [elemIdx: ElementIndex]: Color | undefined,
            },
        },
    }

    /** Provider for custom structure property `SequenceColorProperty` */
    export const Provider: CustomStructureProperty.Provider<Params, Data> = CustomStructureProperty.createProvider({
        label: 'Sequence Color',
        descriptor: CustomPropertyDescriptor<any, any>({
            name: 'sequence-color',
        }),
        type: 'root',
        defaultParams: Params,
        getParams: (data: Structure) => Params,
        isApplicable: (data: Structure) => data.root === data,
        obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<Props>) => {
            const fullProps = { ...PD.getDefaultValues(Params), ...props };
            const items = fullProps.colors.map(t => ({
                selector: t.selector,
                color: t.color,
            } satisfies Data['items'][number]));
            return { value: { items, colorCache: {} } } satisfies CustomProperty.Data<Data>;
        },
    });
}
