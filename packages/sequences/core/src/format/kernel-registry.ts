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

import { KernelAligners, KernelProviders } from "./load/entries.ts";
import { KernelFormats, UniprotService } from "./load/formats.ts";
import { KernelNodeProviders } from "./load/node-providers.ts";
import { KernelLetterMaps } from "./palettes/letter-maps.ts";
import { KernelColorSchemes } from "./palettes/schemes.ts";
import type { SvsRegistryEntry } from "./registry.ts";
import { KernelTables } from "./tables/kernel.ts";

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
