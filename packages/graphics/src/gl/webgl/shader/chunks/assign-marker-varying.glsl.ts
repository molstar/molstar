/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const assign_marker_varying = `
#if defined(dNeedsMarker)
    #if defined(dInstanceExpansion)
        // an expanded object's: its template particle's (every copy's) or-ed with its copy's
        vMarker = float(uint(read_tMarker(granularityIndex(aInstance, group)).r * 255.0 + 0.5) | expansionCopyMarker) / 255.0;
    #else
        vMarker = read_tMarker(granularityIndex(aInstance, group)).r;
    #endif
#endif
`;
