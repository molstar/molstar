/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const clip_pixel = `
#if defined(dClipVariant_pixel) && defined(dClipObjects)
    if (clipTest(vModelPosition / uModelScale))
        discard;
#endif
`;
