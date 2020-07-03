/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { Volume } from '../../mol-model/volume/volume';

export { PropertyProvider };

interface PropertyProvider<T> {
    readonly descriptor: CustomPropertyDescriptor
    get(volume: Volume): T | undefined
    set(volume: Volume, value: T): void
}

namespace PropertyProvider {
    export function create<T>(descriptor: CustomPropertyDescriptor): PropertyProvider<T> {
        const { name } = descriptor;

        return {
            descriptor,
            get(volume: Volume): T | undefined {
                return volume._propertyData[name];
            },
            set(volume: Volume, value: T) {
                volume.customProperties.add(descriptor);
                volume._propertyData[name] = value;
            }
        };
    }
}

//

export { RecommendedIsoValue };

type RecommendedIsoValue = Volume.IsoValue

namespace RecommendedIsoValue {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'recommended_iso_value',
    };

    export const Provider = PropertyProvider.create<RecommendedIsoValue>(Descriptor);
}