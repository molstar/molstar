/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { is_iOS } from '@molstar/core/util/browser';

export const PluginFeatureDetection = {
  get defaultTransparency(): 'blended' | 'wboit' | 'dpoit' {
    return is_iOS() ? 'blended' : 'wboit';
  },
};
