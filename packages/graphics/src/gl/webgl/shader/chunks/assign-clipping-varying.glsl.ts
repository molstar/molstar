/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const assign_clipping_varying = `
#if defined(dClipObjects) && defined(dClipping)
    vClipping = read_tClipping(granularityIndex(aInstance, group)).r;
#endif
`;
