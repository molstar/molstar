/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getMolComponentContext } from './context';
import './index.html';
import './elements/snapshot-markdown';
import './elements/viewer';
import { buildStory } from './kinase-story';
import '../../mol-plugin-ui/skin/light.scss';
import './styles.scss';
import { download } from '../../mol-util/download';

export class MolComponents {
    getContext(name?: string) {
        return getMolComponentContext({ name });
    }
}

(window as any).mc = new MolComponents();
(window as any).buildStory = buildStory;
(window as any).molStarDownload = download;