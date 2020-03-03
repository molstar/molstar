/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor, Model } from '../../../mol-model/structure';
import { ModelFormat } from '../format';

class FormatRegistry<T> {
    map = new Map<ModelFormat['kind'], (model: Model) => T | undefined>()

    add(kind: ModelFormat['kind'], obtain: (model: Model) => T | undefined) {
        this.map.set(kind, obtain)
    }

    remove(kind: ModelFormat['kind']) {
        this.map.delete(kind)
    }

    get(kind: ModelFormat['kind']) {
        return this.map.get(kind)
    }
}

export { FormatPropertyProvider as FormatPropertyProvider }

interface FormatPropertyProvider<T> {
    readonly descriptor: CustomPropertyDescriptor
    readonly formatRegistry: FormatRegistry<T>
    get(model: Model): T | undefined
    set(model: Model, value: T): void
}

namespace FormatPropertyProvider {
    export function create<T>(descriptor: CustomPropertyDescriptor): FormatPropertyProvider<T> {
        const { name } = descriptor
        const formatRegistry = new FormatRegistry<T>()

        return {
            descriptor,
            formatRegistry,
            get(model: Model): T | undefined {
                if (model._staticPropertyData[name]) return model._staticPropertyData[name]
                if (model.customProperties.has(descriptor)) return

                const obtain = formatRegistry.get(model.sourceData.kind)
                if (!obtain) return

                model._staticPropertyData[name] = obtain(model)
                model.customProperties.add(descriptor)
                return model._staticPropertyData[name]
            },
            set(model: Model, value: T) {
                model._staticPropertyData[name] = value
            }
        }
    }
}