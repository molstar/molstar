/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const assign_marker_varying = `
#if defined(dNeedsMarker)
    vMarker = read_tMarker(granularityIndex(aInstance, group)).r;
#endif
`;
