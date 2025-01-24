/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getMolComponentContext } from './context';
import './index.html';
import './snapshot-markdown';
import './viewer';
require('../../mol-plugin-ui/skin/light.scss');

export class MolComponents {
    getContext(name?: string) {
        return getMolComponentContext({ name });
    }
}

(window as any).mc = new MolComponents();