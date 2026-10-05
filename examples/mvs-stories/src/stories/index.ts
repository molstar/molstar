/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { buildStory as kinase } from './kinase.js';
import { buildStory as tbp } from './tbp.js';
import { buildStory as animation } from './animation.js';
import { buildStory as audio } from './audio.js';
import { buildStory as motm1 } from './motm1.js';
import { buildStory as molql } from './molql.js';

export const Stories = [
  { id: 'kinase', name: 'BCR-ABL: A Kinase Out of Control', buildStory: kinase },
  { id: 'tata', name: 'TATA-Binding Protein and its Role in Transcription Initiation ', buildStory: tbp },
  { id: 'motm1', name: 'RCSB PDB Molecule of the Month #1', buildStory: motm1 },
  { id: 'animation-example', name: 'Molecular Animation Example', buildStory: animation },
  { id: 'audio-example', name: 'Audio Playback Example', buildStory: audio },
  { id: 'molql', name: 'Experimental MolQL Selectors', buildStory: molql },
] as const;
