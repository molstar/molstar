/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor, Model, Structure } from '../../mol-model/structure';
import { OrderedMap } from 'immutable';
import { ParamDefinition } from '../../mol-util/param-definition';
import { Task } from '../../mol-task';

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
        return ParamDefinition.MultiSelect(selected, options);
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