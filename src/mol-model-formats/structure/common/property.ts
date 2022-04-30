/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from '../../../mol-model/structure';
import { ModelFormat } from '../../format';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';

class FormatRegistry<T> {
    private map = new Map<ModelFormat['kind'], (model: Model) => T | undefined>();
    private applicable = new Map<ModelFormat['kind'], (model: Model) => boolean>();

    add(kind: ModelFormat['kind'], obtain: (model: Model) => T | undefined, applicable?: (model: Model) => boolean) {
        this.map.set(kind, obtain);
        if (applicable) this.applicable.set(kind, applicable);
    }

    remove(kind: ModelFormat['kind']) {
        this.map.delete(kind);
        this.applicable.delete(kind);
    }

    get(kind: ModelFormat['kind']) {
        return this.map.get(kind);
    }

    isApplicable(model: Model) {
        if (!this.map.has(model.sourceData.kind)) return false;
        const isApplicable = this.applicable.get(model.sourceData.kind);
        return isApplicable ? isApplicable(model) : true;
    }
}

export { FormatPropertyProvider };

interface FormatPropertyProvider<T> {
    readonly descriptor: CustomPropertyDescriptor
    readonly formatRegistry: FormatRegistry<T>
    isApplicable(model: Model): boolean
    get(model: Model): T | undefined
    set(model: Model, value: T): void
    delete(model: Model): void
}

namespace FormatPropertyProvider {
    export function create<T>(descriptor: CustomPropertyDescriptor, options?: { asDynamic?: boolean }): FormatPropertyProvider<T> {
        const { name } = descriptor;
        const formatRegistry = new FormatRegistry<T>();

        return {
            descriptor,
            formatRegistry,
            isApplicable(model: Model) {
                return formatRegistry.isApplicable(model);
            },
            get(model: Model): T | undefined {
                const store = options?.asDynamic ? model._dynamicPropertyData : model._staticPropertyData;

                if (store[name]) return store[name];
                if (model.customProperties.has(descriptor)) return;

                const obtain = formatRegistry.get(model.sourceData.kind);
                if (!obtain) return;

                store[name] = obtain(model);
                model.customProperties.add(descriptor);
                return store[name];
            },
            set(model: Model, value: T) {
                if (options?.asDynamic) {
                    model._dynamicPropertyData[name] = value;
                } else {
                    model._staticPropertyData[name] = value;
                }
            },
            delete(model: Model) {
                if (options?.asDynamic) {
                    delete model._dynamicPropertyData[name];
                } else {
                    delete model._staticPropertyData[name];
                }
            }
        };
    }
}