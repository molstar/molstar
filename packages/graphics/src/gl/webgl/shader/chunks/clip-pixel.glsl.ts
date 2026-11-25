/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const clip_pixel = `
#if defined(dClipVariant_pixel) && defined(dClipObjects)
    // the ghost pass (\`uClipGhost\`) keeps what the clip objects remove
    if (clipTest(vModelPosition / uModelScale) != uClipGhost)
        discard;
#endif
`;
