/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every built-in of seqspec:format#the-document (seqspec:format#extensions): layer types, themes and the `svs.*` codes,
// and the formats, the `uniprot` service, the kernel providers, the node providers (`nj`, `upgma`) and the aligners of
// load/.

import { KernelAligners, KernelProviders } from './load/entries.js';
import { KernelFormats, UniprotService } from './load/formats.js';
import { KernelNodeProviders } from './load/node-providers.js';
import { KernelLetterMaps } from './palettes/letter-maps.js';
import { KernelColorSchemes } from './palettes/schemes.js';
import type { SvsRegistryEntry } from './registry.js';
import { KernelTables } from './tables/kernel.js';

export const KernelRegistry: SvsRegistryEntry = {
  ...KernelTables,
  formats: KernelFormats,
  services: [UniprotService],
  providers: KernelProviders,
  nodeProviders: KernelNodeProviders,
  aligners: KernelAligners,
  letterMaps: KernelLetterMaps,
  colorSchemes: KernelColorSchemes,
};
