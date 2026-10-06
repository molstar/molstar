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
  ['relion_star_particles', RelionStarParticlesProvider] as const,
  ['dynamo_tbl_particles', DynamoTblParticlesProvider] as const,
  ['cryoet_ndjson_particles', CryoEtDataPortalNdjsonParticlesProvider] as const,
  ['artiatomi_em_particles', ArtiatomiEmParticlesProvider] as const,
  ['mmcif_particles', MmcifParticlesProvider] as const,
  ['simularium_particles', SimulariumParticlesProvider] as const,
] as const;

export type BuiltInParticlesFormat = (typeof BuiltInParticlesFormats)[number][0];
