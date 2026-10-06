/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { EmptyPreset } from './empty.js';
import { AutoPreset } from './auto.js';
import { AtomicDetailPreset } from './atomic-detail.js';
import { PolymerCartoonPreset } from './polymer-cartoon.js';
import { PolymerAndLigandPreset } from './polymer-and-ligand.js';
import { ProteinAndNucleicPreset } from './protein-and-nucleic.js';
import { CoarseSurfacePreset } from './coarse-surface.js';
import { IllustrativePreset } from './illustrative.js';
import { MolecularSurfacePreset } from './molecular-surface.js';
import { AutoLodPreset } from './auto-lod.js';
import { MesoscalePreset } from './mesoscale.js';

export const PresetStructureRepresentations = {
  empty: EmptyPreset,
  auto: AutoPreset,
  'atomic-detail': AtomicDetailPreset,
  'polymer-cartoon': PolymerCartoonPreset,
  'polymer-and-ligand': PolymerAndLigandPreset,
  'protein-and-nucleic': ProteinAndNucleicPreset,
  'coarse-surface': CoarseSurfacePreset,
  illustrative: IllustrativePreset,
  'molecular-surface': MolecularSurfacePreset,
  'auto-lod': AutoLodPreset,
  mesoscale: MesoscalePreset,
};
export type PresetStructureRepresentations = typeof PresetStructureRepresentations;

export type BuiltInStructureRepresentationPresetId =
  PresetStructureRepresentations[keyof PresetStructureRepresentations]['id'];
export type BuiltInStructureRepresentationPresetAlias = NonNullable<
  PresetStructureRepresentations[keyof PresetStructureRepresentations]['alias']
>;
