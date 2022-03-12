/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider } from './provider';

export const CoordinatesFormatCategory = 'Coordinates';

export { DcdProvider };
const DcdProvider = DataFormatProvider({
    label: 'DCD',
    description: 'DCD',
    category: CoordinatesFormatCategory,
    binaryExtensions: ['dcd'],
    parse: (plugin, data) => {
        const coordinates = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Model.CoordinatesFromDcd);

        return coordinates.commit();
    }
});
type DcdProvider = typeof DcdProvider;

export { XtcProvider };
const XtcProvider = DataFormatProvider({
    label: 'XTC',
    description: 'XTC',
    category: CoordinatesFormatCategory,
    binaryExtensions: ['xtc'],
    parse: (plugin, data) => {
        const coordinates = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Model.CoordinatesFromXtc);

        return coordinates.commit();
    }
});
type XtcProvider = typeof XtcProvider;

export { TrrProvider };
const TrrProvider = DataFormatProvider({
    label: 'TRR',
    description: 'TRR',
    category: CoordinatesFormatCategory,
    binaryExtensions: ['trr'],
    parse: (plugin, data) => {
        const coordinates = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Model.CoordinatesFromTrr);

        return coordinates.commit();
    }
});
type TrrProvider = typeof TrrProvider;

export type CoordinatesProvider = DcdProvider | XtcProvider | TrrProvider;

export const BuiltInCoordinatesFormats = [
    ['dcd', DcdProvider] as const,
    ['xtc', XtcProvider] as const,
    ['trr', TrrProvider] as const,
] as const;

export type BuiltInCoordinatesFormat = (typeof BuiltInCoordinatesFormats)[number][0]
