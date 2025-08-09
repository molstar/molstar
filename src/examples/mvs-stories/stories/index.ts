/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { buildStory as kinase } from './kinase';
import { buildStory as tbp } from './tbp';
import { buildStory as animation } from './animation'

export const Stories = [
    { id: 'animation', name: 'Molecular Animation', buildStory: animation },
    { id: 'kinase', name: 'BCR-ABL: A Kinase Out of Control', buildStory: kinase },
    { id: 'tata', name: 'TATA-Binding Protein and its Role in Transcription Initiation ', buildStory: tbp },
] as const;