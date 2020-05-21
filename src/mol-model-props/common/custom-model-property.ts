/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueBox } from '../../mol-util';
import { CustomProperty } from './custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { stringToWords } from '../../mol-util/string';

export { CustomModelProperty };

namespace CustomModelProperty {
    export interface Provider<Params extends PD.Params, Value> extends CustomProperty.Provider<Model, Params, Value> { }

    export interface ProviderBuilder<Params extends PD.Params, Value> {
        readonly label: string
        readonly descriptor: CustomPropertyDescriptor
        readonly isHidden?: boolean
        readonly defaultParams: Params
        readonly getParams: (data: Model) => Params
        readonly isApplicable: (data: Model) => boolean
        readonly obtain: (ctx: CustomProperty.Context, data: Model, props: PD.Values<Params>) => Promise<CustomProperty.Data<Value>>
        readonly type: 'static' | 'dynamic'
    }

    export function createProvider<Params extends PD.Params, Value>(builder: ProviderBuilder<Params, Value>): CustomProperty.Provider<Model, Params, Value> {
        const descriptorName = builder.descriptor.name;
        const propertyDataName = builder.type === 'static' ? '_staticPropertyData' : '_dynamicPropertyData';

        const get = (data: Model) => {
            if (!(descriptorName in data[propertyDataName])) {
                (data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>) = {
                    props: { ...PD.getDefaultValues(builder.getParams(data)) },
                    data: ValueBox.create(undefined)
                };
            }
            return data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>;
        };
        const set = (data: Model, props: PD.Values<Params>, value: Value | undefined) => {
            const property = get(data);
            (data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>) = {
                props,
                data: ValueBox.withValue(property.data, value)
            };
        };

        return {
            label: builder.label,
            descriptor: builder.descriptor,
            isHidden: builder.isHidden,
            getParams: (data: Model) => {
                const params = PD.clone(builder.getParams(data));
                PD.setDefaultValues(params, get(data).props);
                return params;
            },
            defaultParams: builder.defaultParams,
            isApplicable: builder.isApplicable,
            attach: async (ctx: CustomProperty.Context, data: Model, props: Partial<PD.Values<Params>> = {}, addRef) => {
                if (addRef) data.customProperties.reference(builder.descriptor, true);
                const property = get(data);
                const p = PD.merge(builder.defaultParams, property.props, props);
                if (property.data.value && PD.areEqual(builder.defaultParams, property.props, p)) return;
                const { value, assets } = await builder.obtain(ctx, data, p);
                data.customProperties.add(builder.descriptor);
                data.customProperties.assets(builder.descriptor, assets);
                set(data, p, value);
            },
            ref: (data: Model, add: boolean) => data.customProperties.reference(builder.descriptor, add),
            get: (data: Model) => get(data)?.data,
            set: (data: Model, props: Partial<PD.Values<Params>> = {}) => {
                const property = get(data);
                const p = PD.merge(builder.defaultParams, property.props, props);
                if (!PD.areEqual(builder.defaultParams, property.props, p)) {
                    // this invalidates property.value
                    set(data, p, undefined);
                    // dispose of assets
                    data.customProperties.assets(builder.descriptor);
                }
            },
            props: (data: Model) => get(data).props,
        };
    }

    export function createSimple<T>(name: string, type: 'static' | 'dynamic', defaultValue?: T) {
        const defaultParams = { value: PD.Value(defaultValue, { isHidden: true }) };
        return createProvider({
            label: stringToWords(name),
            descriptor: CustomPropertyDescriptor({ name }),
            isHidden: true,
            type,
            defaultParams,
            getParams: () => ({ value: PD.Value(defaultValue, { isHidden: true }) }),
            isApplicable: () => true,
            obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<PD.Values<typeof defaultParams>>) => {
                return { value: props.value ?? defaultValue };
            }
        });
    }
}