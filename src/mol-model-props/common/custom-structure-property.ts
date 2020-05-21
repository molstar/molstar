/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueBox } from '../../mol-util';
import { CustomProperty } from './custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { stringToWords } from '../../mol-util/string';

export { CustomStructureProperty };

namespace CustomStructureProperty {
    export interface Provider<Params extends PD.Params, Value> extends CustomProperty.Provider<Structure, Params, Value> { }

    export interface ProviderBuilder<Params extends PD.Params, Value> {
        readonly label: string
        readonly descriptor: CustomPropertyDescriptor
        readonly isHidden?: boolean
        readonly defaultParams: Params
        readonly getParams: (data: Structure) => Params
        readonly isApplicable: (data: Structure) => boolean
        readonly obtain: (ctx: CustomProperty.Context, data: Structure, props: PD.Values<Params>) => Promise<CustomProperty.Data<Value>>
        readonly type: 'root' | 'local'
    }

    export function createProvider<Params extends PD.Params, Value>(builder: ProviderBuilder<Params, Value>): CustomProperty.Provider<Structure, Params, Value> {
        const descriptorName = builder.descriptor.name;
        const propertyDataName = builder.type === 'root' ? 'inheritedPropertyData' : 'currentPropertyData';

        const get = (data: Structure) => {
            if (!(descriptorName in data[propertyDataName])) {
                (data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>) = {
                    props: { ...PD.getDefaultValues(builder.getParams(data)) },
                    data: ValueBox.create(undefined)
                };
            }
            return data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>;
        };
        const set = (data: Structure, props: PD.Values<Params>, value: Value | undefined) => {
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
            getParams: (data: Structure) => {
                const params = PD.clone(builder.getParams(data));
                PD.setDefaultValues(params, get(data).props);
                return params;
            },
            defaultParams: builder.defaultParams,
            isApplicable: builder.isApplicable,
            attach: async (ctx: CustomProperty.Context, data: Structure, props: Partial<PD.Values<Params>> = {}, addRef) => {
                if (addRef) data.customPropertyDescriptors.reference(builder.descriptor, true);
                if (builder.type === 'root') data = data.root;
                const rootProps = get(data.root).props;
                const property = get(data);
                const p = PD.merge(builder.defaultParams, rootProps, props);
                if (property.data.value && PD.areEqual(builder.defaultParams, property.props, p)) return;
                const { value, assets } = await builder.obtain(ctx, data, p);
                data.customPropertyDescriptors.add(builder.descriptor);
                data.customPropertyDescriptors.assets(builder.descriptor, assets);
                set(data, p, value);
            },
            ref: (data: Structure, add: boolean) => data.customPropertyDescriptors.reference(builder.descriptor, add),
            get: (data: Structure) => get(data).data,
            set: (data: Structure, props: Partial<PD.Values<Params>> = {}, value?: Value) => {
                if (builder.type === 'root') data = data.root;
                const property = get(data);
                const p = PD.merge(builder.defaultParams, property.props, props);
                if (!PD.areEqual(builder.defaultParams, property.props, p)) {
                    // this invalidates property.value
                    set(data, p, value);
                    // dispose of assets
                    data.customPropertyDescriptors.assets(builder.descriptor);
                }
            },
            props: (data: Structure) => get(data).props,
        };
    }

    export function createSimple<T>(name: string, type: 'root' | 'local', defaultValue?: T) {
        const defaultParams = { value: PD.Value(defaultValue, { isHidden: true }) };
        return createProvider({
            label: stringToWords(name),
            descriptor: CustomPropertyDescriptor({ name }),
            isHidden: true,
            type,
            defaultParams,
            getParams: () => ({ value: PD.Value(defaultValue, { isHidden: true }) }),
            isApplicable: () => true,
            obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<PD.Values<typeof defaultParams>>) => {
                return { ...PD.getDefaultValues(defaultParams), ...props };
            }
        });
    }
}