/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../mol-model/structure';
import { ModelFormat } from '../../format';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';

class FormatRegistry<T> {
    private map = new Map<ModelFormat['kind'], (model: Model) => T | undefined>()
    private applicable = new Map<ModelFormat['kind'], (model: Model) => boolean>()

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
        const isApplicable = this.applicable.get(model.sourceData.kind);
        return isApplicable ? isApplicable(model) : true;
    }
}

export { FormatPropertyProvider as FormatPropertyProvider };

interface FormatPropertyProvider<T> {
    readonly descriptor: CustomPropertyDescriptor
    readonly formatRegistry: FormatRegistry<T>
    isApplicable(model: Model): boolean
    get(model: Model): T | undefined
    set(model: Model, value: T): void
}

namespace FormatPropertyProvider {
    export function create<T>(descriptor: CustomPropertyDescriptor): FormatPropertyProvider<T> {
        const { name } = descriptor;
        const formatRegistry = new FormatRegistry<T>();

        return {
            descriptor,
            formatRegistry,
            isApplicable(model: Model) {
                return formatRegistry.isApplicable(model);
            },
            get(model: Model): T | undefined {
                if (model._staticPropertyData[name]) return model._staticPropertyData[name];
                if (model.customProperties.has(descriptor)) return;

                const obtain = formatRegistry.get(model.sourceData.kind);
                if (!obtain) return;

                model._staticPropertyData[name] = obtain(model);
                model.customProperties.add(descriptor);
                return model._staticPropertyData[name];
            },
            set(model: Model, value: T) {
                model._staticPropertyData[name] = value;
            }
        };
    }
}