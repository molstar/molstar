/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { namedCatalog } from '@molstar/graphics/util/named-catalog';
import { AtomIdColorThemeProvider } from './atom-id.js';
import { CarbohydrateSymbolColorThemeProvider } from './carbohydrate-symbol.js';
import { CartoonColorThemeProvider } from './cartoon.js';
import { ChainIdColorThemeProvider } from './chain-id.js';
import { ElementIndexColorThemeProvider } from './element-index.js';
import { ElementSymbolColorThemeProvider } from './element-symbol.js';
import { EntityIdColorThemeProvider } from './entity-id.js';
import { EntitySourceColorThemeProvider } from './entity-source.js';
import { FormalChargeColorThemeProvider } from './formal-charge.js';
import { HydrophobicityColorThemeProvider } from './hydrophobicity.js';
import { IllustrativeColorThemeProvider } from './illustrative.js';
import { ModelIndexColorThemeProvider } from './model-index.js';
import { MoleculeTypeColorThemeProvider } from './molecule-type.js';
import { OccupancyColorThemeProvider } from './occupancy.js';
import { OperatorHklColorThemeProvider } from './operator-hkl.js';
import { OperatorNameColorThemeProvider } from './operator-name.js';
import { PartialChargeColorThemeProvider } from './partial-charge.js';
import { ParticleAttributeColorThemeProvider } from './particle-attribute.js';
import { ParticleCompartmentColorThemeProvider } from './particle-compartment.js';
import { ParticleEntityColorThemeProvider } from './particle-entity.js';
import { ParticleHierarchyColorThemeProvider } from './particle-hierarchy.js';
import { ParticleIndexColorThemeProvider } from './particle-index.js';
import { PolymerIdColorThemeProvider } from './polymer-id.js';
import { PolymerIndexColorThemeProvider } from './polymer-index.js';
import { ResidueChargeColorThemeProvider } from './residue-charge.js';
import { ResidueNameColorThemeProvider } from './residue-name.js';
import { SecondaryStructureColorThemeProvider } from './secondary-structure.js';
import { SequenceIdColorThemeProvider } from './sequence-id.js';
import { ShapeGroupColorThemeProvider } from './shape-group.js';
import { StructureIndexColorThemeProvider } from './structure-index.js';
import { TrajectoryIndexColorThemeProvider } from './trajectory-index.js';
import { UncertaintyColorThemeProvider } from './uncertainty.js';
import { UnitIndexColorThemeProvider } from './unit-index.js';
import { UniformColorThemeProvider } from './uniform.js';
import { VolumeInstanceColorThemeProvider } from './volume-instance.js';
import { VolumeSegmentColorThemeProvider } from './volume-segment.js';
import { VolumeValueColorThemeProvider } from './volume-value.js';

export const BuiltInColorThemes = namedCatalog({
  'atom-id': AtomIdColorThemeProvider,
  'carbohydrate-symbol': CarbohydrateSymbolColorThemeProvider,
  cartoon: CartoonColorThemeProvider,
  'chain-id': ChainIdColorThemeProvider,
  'element-index': ElementIndexColorThemeProvider,
  'element-symbol': ElementSymbolColorThemeProvider,
  'entity-id': EntityIdColorThemeProvider,
  'entity-source': EntitySourceColorThemeProvider,
  'formal-charge': FormalChargeColorThemeProvider,
  hydrophobicity: HydrophobicityColorThemeProvider,
  illustrative: IllustrativeColorThemeProvider,
  'model-index': ModelIndexColorThemeProvider,
  'molecule-type': MoleculeTypeColorThemeProvider,
  occupancy: OccupancyColorThemeProvider,
  'operator-hkl': OperatorHklColorThemeProvider,
  'operator-name': OperatorNameColorThemeProvider,
  'partial-charge': PartialChargeColorThemeProvider,
  'particle-attribute': ParticleAttributeColorThemeProvider,
  'particle-compartment': ParticleCompartmentColorThemeProvider,
  'particle-entity': ParticleEntityColorThemeProvider,
  'particle-hierarchy': ParticleHierarchyColorThemeProvider,
  'particle-index': ParticleIndexColorThemeProvider,
  'polymer-id': PolymerIdColorThemeProvider,
  'polymer-index': PolymerIndexColorThemeProvider,
  'residue-charge': ResidueChargeColorThemeProvider,
  'residue-name': ResidueNameColorThemeProvider,
  'secondary-structure': SecondaryStructureColorThemeProvider,
  'sequence-id': SequenceIdColorThemeProvider,
  'shape-group': ShapeGroupColorThemeProvider,
  'structure-index': StructureIndexColorThemeProvider,
  'trajectory-index': TrajectoryIndexColorThemeProvider,
  uncertainty: UncertaintyColorThemeProvider,
  'unit-index': UnitIndexColorThemeProvider,
  uniform: UniformColorThemeProvider,
  'volume-instance': VolumeInstanceColorThemeProvider,
  'volume-segment': VolumeSegmentColorThemeProvider,
  'volume-value': VolumeValueColorThemeProvider,
});
