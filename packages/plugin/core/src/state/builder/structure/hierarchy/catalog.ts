/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DefaultHierarchyPreset } from './default.js';
import { AllModelsHierarchyPreset } from './all-models.js';
import { UnitcellHierarchyPreset } from './unitcell.js';
import { SupercellHierarchyPreset } from './supercell.js';
import { CrystalContactsHierarchyPreset } from './crystal-contacts.js';

export const PresetTrajectoryHierarchy = {
  default: DefaultHierarchyPreset,
  'all-models': AllModelsHierarchyPreset,
  unitcell: UnitcellHierarchyPreset,
  supercell: SupercellHierarchyPreset,
  crystalContacts: CrystalContactsHierarchyPreset,
};
export type PresetTrajectoryHierarchy = typeof PresetTrajectoryHierarchy;

export type BuiltInTrajectoryHierarchyPresetId = PresetTrajectoryHierarchy[keyof PresetTrajectoryHierarchy]['id'];
export type BuiltInTrajectoryHierarchyPresetAlias = NonNullable<
  PresetTrajectoryHierarchy[keyof PresetTrajectoryHierarchy]['alias']
>;
