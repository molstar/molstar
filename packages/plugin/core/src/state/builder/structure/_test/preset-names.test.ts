/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type BuiltInTrajectoryHierarchyPresetAlias,
  type BuiltInTrajectoryHierarchyPresetId,
  PresetTrajectoryHierarchy,
} from '../hierarchy-presets/catalog.js';
import {
  type BuiltInStructureRepresentationPresetAlias,
  type BuiltInStructureRepresentationPresetId,
  PresetStructureRepresentations,
} from '../representation-presets/catalog.js';

type Equal<A, B> = (<T>() => T extends A ? 1 : 2) extends <T>() => T extends B ? 1 : 2 ? true : false;
function assertType<_T extends true>() {}

assertType<
  Equal<
    BuiltInTrajectoryHierarchyPresetId,
    | 'preset-trajectory-default'
    | 'preset-trajectory-all-models'
    | 'preset-trajectory-unitcell'
    | 'preset-trajectory-supercell'
    | 'preset-trajectory-crystal-contacts'
  >
>();
assertType<
  Equal<BuiltInTrajectoryHierarchyPresetAlias, 'default' | 'all-models' | 'unitcell' | 'supercell' | 'crystalContacts'>
>();
assertType<Equal<BuiltInTrajectoryHierarchyPresetAlias, keyof PresetTrajectoryHierarchy>>();

assertType<
  Equal<
    BuiltInStructureRepresentationPresetId,
    | 'preset-structure-representation-empty'
    | 'preset-structure-representation-auto'
    | 'preset-structure-representation-atomic-detail'
    | 'preset-structure-representation-polymer-cartoon'
    | 'preset-structure-representation-polymer-and-ligand'
    | 'preset-structure-representation-protein-and-nucleic'
    | 'preset-structure-representation-coarse-surface'
    | 'preset-structure-representation-illustrative'
    | 'preset-structure-representation-molecular-surface'
    | 'preset-structure-representation-auto-lod'
    | 'preset-structure-representation-mesoscale'
  >
>();
assertType<
  Equal<
    BuiltInStructureRepresentationPresetAlias,
    | 'empty'
    | 'auto'
    | 'atomic-detail'
    | 'polymer-cartoon'
    | 'polymer-and-ligand'
    | 'protein-and-nucleic'
    | 'coarse-surface'
    | 'illustrative'
    | 'molecular-surface'
    | 'auto-lod'
    | 'mesoscale'
  >
>();
assertType<Equal<BuiltInStructureRepresentationPresetAlias, keyof PresetStructureRepresentations>>();

describe('built-in preset names', () => {
  it('sets the alias of every hierarchy preset to its catalog key', () => {
    for (const [key, preset] of Object.entries(PresetTrajectoryHierarchy)) {
      expect(preset.alias).toBe(key);
    }
  });

  it('sets the alias of every structure representation preset to its catalog key', () => {
    for (const [key, preset] of Object.entries(PresetStructureRepresentations)) {
      expect(preset.alias).toBe(key);
      expect(preset.id).toBe(`preset-structure-representation-${key}`);
    }
  });

  it('keeps the catalog order', () => {
    expect(Object.keys(PresetTrajectoryHierarchy)).toEqual([
      'default',
      'all-models',
      'unitcell',
      'supercell',
      'crystalContacts',
    ]);
    expect(Object.keys(PresetStructureRepresentations)).toEqual([
      'empty',
      'auto',
      'atomic-detail',
      'polymer-cartoon',
      'polymer-and-ligand',
      'protein-and-nucleic',
      'coarse-surface',
      'illustrative',
      'molecular-surface',
      'auto-lod',
      'mesoscale',
    ]);
  });
});
