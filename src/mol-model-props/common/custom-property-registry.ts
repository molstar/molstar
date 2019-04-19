/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CustomPropertyDescriptor, Model } from 'mol-model/structure';
import { OrderedMap } from 'immutable';
import { ParamDefinition } from 'mol-util/param-definition';
import { Task } from 'mol-task';

export { CustomPropertyRegistry }

class CustomPropertyRegistry {
    private providers = OrderedMap<string, CustomPropertyRegistry.Provider>().asMutable();

    getSelect(model: Model) {
        const values = this.providers.values();
        const options: [string, string][] = [], selected: string[] = [];
        while (true) {
            const v = values.next();
            if (v.done) break;
            if (!v.value.attachableTo(model)) continue;
            options.push(v.value.option);
            if (v.value.defaultSelected) selected.push(v.value.option[0]);
        }
        return ParamDefinition.MultiSelect(selected, options);
    }

    getDefault(model: Model) {
        const values = this.providers.values();
        const selected: string[] = [];
        while (true) {
            const v = values.next();
            if (v.done) break;
            if (!v.value.attachableTo(model)) continue;
            if (v.value.defaultSelected) selected.push(v.value.option[0]);
        }
        return selected;
    }

    get(name: string) {
        const prop = this.providers.get(name);
        if (!prop) throw new Error(`Custom prop '${name}' is not registered.`);
        return this.providers.get(name);
    }

    register(provider: CustomPropertyRegistry.Provider) {
        this.providers.set(provider.descriptor.name, provider);
    }

    unregister(name: string) {
        this.providers.delete(name);
    }
}

namespace CustomPropertyRegistry {
    export interface Provider {
        option: [string, string],
        defaultSelected: boolean,
        descriptor: CustomPropertyDescriptor<any, any>,
        attachableTo: (model: Model) => boolean,
        attach: (model: Model) => Task<boolean>
    }
}