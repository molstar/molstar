/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome chunk (_spec/README.md "Resolve", Assemblies, Loaded when used): imported by `loadGenome` only.

import { assemblyEntry } from './assemblies.js';
import { assemblyAnnotation } from './genome-data.js';
import { assemblyPairs, checkAssemblyData } from './genome-model.js';
import { assemblyLengths, listedSystems, resolveAssemblies } from './genome-resolve.js';

export const Genome = {
  assemblyAnnotation,
  assemblyEntry,
  assemblyLengths,
  assemblyPairs,
  checkAssemblyData,
  listedSystems,
  resolveAssemblies,
};
