/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { buildStory as kinase } from './kinase';
import { buildStory as tbp } from './tbp';
import { buildStory as animation } from './animation';
import { buildStory as trajectory } from './trajectory';

export const Stories = [
    { id: 'kinase', name: 'BCR-ABL: A Kinase Out of Control', buildStory: kinase },
    { id: 'tata', name: 'TATA-Binding Protein and its Role in Transcription Initiation ', buildStory: tbp },
    { id: 'animation', name: 'Molecular Animation', buildStory: animation },
    // Temp
    { id: 'traj', name: 'Trajectory', buildStory: trajectory },
] as const;