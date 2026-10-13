/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RelionStarParticlesProvider } from './star.js';
import { DynamoTblParticlesProvider } from './tbl.js';
import { CryoEtDataPortalNdjsonParticlesProvider } from './ndjson.js';
import { ArtiatomiEmParticlesProvider } from './em.js';
import { MmcifParticlesProvider } from './mmcif-assembly.js';
import { SimulariumParticlesProvider } from './simularium.js';

export const BuiltInParticlesFormats = [
  RelionStarParticlesProvider,
  DynamoTblParticlesProvider,
  CryoEtDataPortalNdjsonParticlesProvider,
  ArtiatomiEmParticlesProvider,
  MmcifParticlesProvider,
  SimulariumParticlesProvider,
] as const;

export type BuiltInParticlesFormat = (typeof BuiltInParticlesFormats)[number]['name'];
