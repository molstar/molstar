/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor, Model, Structure } from '../../mol-model/structure';
import { OrderedMap } from 'immutable';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Task, RuntimeContext } from '../../mol-task';
import { ValueBox } from '../../mol-util';

export { CustomPropertyRegistry }

class CustomPropertyRegistry<T = never> {
    private providers = OrderedMap<string, CustomPropertyRegistry.Provider<T>>().asMutable();

    getSelect(object: T) {
        const values = this.providers.values();
        const options: [string, string][] = [], selected: string[] = [];
        while (true) {
            const v = values.next();
            if (v.done) break;
            if (!v.value.attachableTo(object)) continue;
            options.push(v.value.option);
            if (v.value.defaultSelected) selected.push(v.value.option[0]);
        }
        return PD.MultiSelect(selected, options);
    }

    getDefault(object: T) {
        const values = this.providers.values();
        const selected: string[] = [];
        while (true) {
            const v = values.next();
            if (v.done) break;
            if (!v.value.attachableTo(object)) continue;
            if (v.value.defaultSelected) selected.push(v.value.option[0]);
        }
        return selected;
    }

    get(name: string) {
        const prop = this.providers.get(name);
        if (!prop) throw new Error(`Custom prop '${name}' is not registered.`);
        return this.providers.get(name);
    }

    register(provider: CustomPropertyRegistry.Provider<T>) {
        this.providers.set(provider.descriptor.name, provider);
    }

    unregister(name: string) {
        this.providers.delete(name);
    }
}

namespace CustomPropertyRegistry {
    export interface Provider<T> {
        option: [string, string],
        defaultSelected: boolean,
        descriptor: CustomPropertyDescriptor<any, any>,
        attachableTo: (object: T) => boolean,
        attach: (object: T) => Task<boolean>
    }

    export type ModelProvider = Provider<Model>
    export type StructureProvider = Provider<Structure>
}

export { CustomStructureProperty }

namespace CustomStructureProperty {
    export interface Provider<Params extends PD.Params, Value> {
        label: string
        descriptor: CustomPropertyDescriptor
        getParams: (data: Structure) => Params
        isApplicable: (data: Structure) => boolean
        attach: (data: Structure, props?: Partial<PD.Values<Params>>) => Task<void>
        getValue: (data: Structure) => ValueBox<Value | undefined>
        setProps: (data: Structure, props: PD.Values<Params>) => void
    }

    export interface ProviderBuilder<Params extends PD.Params, Value> {
        label: string
        defaultParams: Params
        getParams: (data: Structure) => Params
        isApplicable: (data: Structure) => boolean
        compute: (ctx: RuntimeContext, data: Structure, props: PD.Values<Params>) => Promise<Value>
        descriptor: CustomPropertyDescriptor
    }

    // TODO currently this always uses .inheritedPropertyData
    export function createProvider<Params extends PD.Params, Value>(builder: ProviderBuilder<Params, Value>): CustomStructureProperty.Provider<Params, Value> {
        const get = (data: Structure) => {
            if (!(builder.descriptor.name in data.inheritedPropertyData)) {
                (data.inheritedPropertyData[builder.descriptor.name] as CustomStructureProperty.Property<PD.Values<Params>, Value>) = {
                    props: { ...PD.getDefaultValues(builder.getParams(data)) },
                    data: ValueBox.create(undefined)
                }
            }
            return data.inheritedPropertyData[builder.descriptor.name] as CustomStructureProperty.Property<PD.Values<Params>, Value>;
        }
        const set = (data: Structure, props: PD.Values<Params>, value: Value | undefined) => {
            const property = get(data);
            (data.inheritedPropertyData[builder.descriptor.name] as CustomStructureProperty.Property<PD.Values<Params>, Value>) = {
                props,
                data: ValueBox.withValue(property.data, value)
            };
        }

        return {
            label: builder.label,
            descriptor: builder.descriptor,
            getParams: builder.getParams,
            isApplicable: builder.isApplicable,
            attach: (data: Structure, props: Partial<PD.Values<Params>> = {}) => Task.create(`Attach ${builder.label}`, async ctx => {
                const property = get(data)
                const p = { ...property.props, ...props }
                if (property.data.value && PD.areEqual(builder.defaultParams, property.props, p)) return
                const value = await builder.compute(ctx, data, p)
                data.customPropertyDescriptors.add(builder.descriptor);
                set(data, p, value);
            }),
            getValue: (data: Structure) => get(data)?.data,
            setProps: (data: Structure, props: Partial<PD.Values<Params>> = {}) => {
                const property = get(data)
                const p = { ...property.props, ...props }
                if (!PD.areEqual(builder.defaultParams, property.props, p)) {
                    // this invalidates property.value
                    set(data, p, undefined)
                }
            },
        }
    }

    export interface Property<P, V> {
        readonly props: P
        readonly data: ValueBox<V | undefined>
    }

    export class Registry {
        private providers = OrderedMap<string, Provider<any, any>>().asMutable();

        /** Get params for all applicable property providers */
        getParams(data: Structure) {
            const values = this.providers.values();
            const params: PD.Params = {};
            while (true) {
                const v = values.next();
                if (v.done) break;
                if (!v.value.isApplicable(data)) continue;
                params[v.value.descriptor.name] = PD.Group(v.value.getParams(data), { label: v.value.label })
            }
            return params
        }

        get(name: string) {
            const prop = this.providers.get(name);
            if (!prop) throw new Error(`Custom property '${name}' is not registered.`);
            return this.providers.get(name);
        }

        register(provider: Provider<any, any>) {
            this.providers.set(provider.descriptor.name, provider);
        }

        unregister(name: string) {
            this.providers.delete(name);
        }
    }
}