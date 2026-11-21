/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome chunk (_spec/README.md "Resolve", Assemblies, Loaded when used): imported by `loadGenome` only.

import { assemblyEntry } from './assemblies.js';
import { assemblyPairs, checkAssemblyData } from './genome-model.js';

export const Genome = { assemblyEntry, assemblyPairs, checkAssemblyData };
