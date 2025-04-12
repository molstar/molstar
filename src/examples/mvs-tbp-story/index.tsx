/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { getMolComponentContext } from '../mvs-kinase-story/context';
import './index.html';
import '../mvs-kinase-story/elements/snapshot-markdown';
import '../mvs-kinase-story/elements/viewer';
import { buildStory } from './tbp-story';
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