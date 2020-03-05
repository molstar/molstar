/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor, Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueBox } from '../../mol-util';
import { CustomProperty } from './custom-property';

export { CustomStructureProperty }

namespace CustomStructureProperty {
    export interface Provider<Params extends PD.Params, Value> extends CustomProperty.Provider<Structure, Params, Value> { }

    export interface ProviderBuilder<Params extends PD.Params, Value> {
        readonly label: string
        readonly descriptor: CustomPropertyDescriptor
        readonly defaultParams: Params
        readonly getParams: (data: Structure) => Params
        readonly isApplicable: (data: Structure) => boolean
        readonly obtain: (ctx: CustomProperty.Context, data: Structure, props: PD.Values<Params>) => Promise<Value>
        readonly type: 'root' | 'local'
    }

    export function createProvider<Params extends PD.Params, Value>(builder: ProviderBuilder<Params, Value>): CustomProperty.Provider<Structure, Params, Value> {
        const descriptorName = builder.descriptor.name
        const propertyDataName = builder.type === 'root' ? 'inheritedPropertyData' : 'currentPropertyData'

        const get = (data: Structure) => {
            if (!(descriptorName in data[propertyDataName])) {
                (data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>) = {
                    props: { ...PD.getDefaultValues(builder.getParams(data)) },
                    data: ValueBox.create(undefined)
                }
            }
            return data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>;
        }
        const set = (data: Structure, props: PD.Values<Params>, value: Value | undefined) => {
            const property = get(data);
            (data[propertyDataName][descriptorName] as CustomProperty.Container<PD.Values<Params>, Value>) = {
                props,
                data: ValueBox.withValue(property.data, value)
            };
        }

        return {
            label: builder.label,
            descriptor: builder.descriptor,
            getParams: builder.getParams,
            defaultParams: builder.defaultParams,
            isApplicable: builder.isApplicable,
            attach: async (ctx: CustomProperty.Context, data: Structure, props: Partial<PD.Values<Params>> = {}) => {
                if (builder.type === 'root') data = data.root
                const property = get(data)
                const p = { ...property.props, ...props }
                if (property.data.value && PD.areEqual(builder.defaultParams, property.props, p)) return
                const value = await builder.obtain(ctx, data, p)
                data.customPropertyDescriptors.add(builder.descriptor);
                set(data, p, value);
            },
            get: (data: Structure) => get(data).data,
            set: (data: Structure, props: Partial<PD.Values<Params>> = {}, value?: Value) => {
                if (builder.type === 'root') data = data.root
                const property = get(data)
                const p = { ...property.props, ...props }
                if (!PD.areEqual(builder.defaultParams, property.props, p)) {
                    // this invalidates property.value
                    set(data, p, value)
                }
            }
        }
    }
}