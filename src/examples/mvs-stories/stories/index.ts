/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { buildStory as kinase } from './kinase';
import { buildStory as tbp } from './tbp';
import { buildStory as animation } from './animation';
import { buildStory as audio } from './audio';
import { buildStory as motm1 } from './motm1';

export const Stories = [
    { id: 'kinase', name: 'BCR-ABL: A Kinase Out of Control', buildStory: kinase },
    { id: 'tata', name: 'TATA-Binding Protein and its Role in Transcription Initiation ', buildStory: tbp },
    { id: 'motm1', name: 'RCSB Molecule of the Month #1', buildStory: motm1 },
    { id: 'animation-example', name: 'Molecular Animation Example', buildStory: animation },
    { id: 'audio-example', name: 'Audio Playback Example', buildStory: audio },
] as const;