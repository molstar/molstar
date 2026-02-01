/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueBox } from '../../mol-util';
import { CustomProperty } from './custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { stringToWords } from '../../mol-util/string';
import { Volume } from '../../mol-model/volume/volume';

export { CustomVolumeProperty };

namespace CustomVolumeProperty {
    export interface Provider<Params extends PD.Params, Value> extends CustomProperty.Provider<Volume, Params, Value> { }

    export interface ProviderBuilder<Params extends PD.Params, Value> {
        readonly label: string
        readonly descriptor: CustomPropertyDescriptor
        readonly isHidden?: boolean
        readonly defaultParams: Params
        readonly getParams: (data: Volume) => Params
        readonly isApplicable: (data: Volume) => boolean
        readonly obtain: (ctx: CustomProperty.Context, data: Volume, props: PD.Values<Params>) => Promise<CustomProperty.Data<Value>>
    }

    export function createProvider<Params extends PD.Params, Value>(builder: ProviderBuilder<Params, Value>): CustomProperty.Provider<Volume, Params, Value> {
        const descriptorName = builder.descriptor.name;

        const get = (data: Volume) => {
            if (!(descriptorName in data._propertyData)) {
                (data._propertyData[descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>) = {
                    props: { ...PD.getDefaultValues(builder.getParams(data)) },
                    data: ValueBox.create(undefined)
                };
            }
            return data._propertyData[descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>;
        };
        const set = (data: Volume, props: PD.Values<Params>, value: Value | undefined) => {
            const property = get(data);
            (data._propertyData[descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>) = {
                props,
                data: ValueBox.withValue(property.data, value)
            };
        };

        return {
            label: builder.label,
            descriptor: builder.descriptor,
            isHidden: builder.isHidden,
            getParams: (data: Volume) => {
                const params = PD.clone(builder.getParams(data));
                PD.setDefaultValues(params, get(data).props);
                return params;
            },
            defaultParams: builder.defaultParams,
            isApplicable: builder.isApplicable,
            attach: async (ctx: CustomProperty.Context, data: Volume, props: Partial<PD.Values<Params>> = {}, addRef) => {
                if (addRef) data.customProperties.reference(builder.descriptor, true);
                const property = get(data);
                const p = PD.merge(builder.defaultParams, property.props, props);
                if (property.data.value && PD.areEqual(builder.defaultParams, property.props, p)) return;
                const { value, assets } = await builder.obtain(ctx, data, p);
                data.customProperties.add(builder.descriptor);
                data.customProperties.assets(builder.descriptor, assets);
                set(data, p, value);
            },
            ref: (data: Volume, add: boolean) => data.customProperties.reference(builder.descriptor, add),
            get: (data: Volume) => get(data)?.data,
            set: (data: Volume, props: Partial<PD.Values<Params>> = {}, value?: Value) => {
                const property = get(data);
                const p = PD.merge(builder.defaultParams, property.props, props);
                if (!PD.areEqual(builder.defaultParams, property.props, p)) {
                    // this invalidates property.value
                    set(data, p, value);
                    // dispose of assets
                    data.customProperties.assets(builder.descriptor);
                }
            },
            props: (data: Volume) => get(data).props,
        };
    }

    export function createSimple<T>(name: string, defaultValue?: T) {
        const defaultParams = { value: PD.Value(defaultValue, { isHidden: true }) };
        return createProvider({
            label: stringToWords(name),
            descriptor: CustomPropertyDescriptor({ name }),
            isHidden: true,
            defaultParams,
            getParams: () => ({ value: PD.Value(defaultValue, { isHidden: true }) }),
            isApplicable: () => true,
            obtain: async (ctx: CustomProperty.Context, data: Volume, props: Partial<PD.Values<typeof defaultParams>>) => {
                return { ...PD.getDefaultValues(defaultParams), ...props };
            }
        });
    }
}