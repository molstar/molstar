/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Russ Taylor <russ@reliasolve.com>
 */

import { ANVILMembraneOrientation } from '@molstar/anvil-extension/behavior';
import { AssemblySymmetry } from '@molstar/assembly-symmetry-extension';
import { Backgrounds } from '@molstar/backgrounds-extension';
import { DebugHelpers } from '@molstar/debug-helpers-extension';
import { DnatcoNtCs } from '@molstar/dnatco-extension';
import { G3DFormat } from '@molstar/g3d-extension/format';
import { GeometryExport } from '@molstar/geo-export-extension';
import {
  MAQualityAssessment,
  MAQualityAssessmentConfig,
} from '@molstar/model-archive-extension/quality-assessment/behavior';
import { ModelExport } from '@molstar/model-export-extension';
import { Mp4Export } from '@molstar/mp4-export-extension';
import { loadMVS } from '@molstar/mvs';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { loadMVSData } from '@molstar/mvs/components/formats';
import { PDBeStructureQualityReport } from '@molstar/pdbe-extension';
import { RCSBValidationReport } from '@molstar/rcsb-extension';
import { SbNcbrPartialCharges, SbNcbrTunnels } from '@molstar/sb-ncbr-extension';
import { wwPDBChemicalComponentDictionary } from '@molstar/wwpdb-extension/ccd/behavior';
import { wwPDBStructConnExtensionFunctions } from '@molstar/wwpdb-extension/struct-conn';
import { ZenodoImport } from '@molstar/zenodo-extension';
import { PluginSpec } from '@molstar/plugin/spec';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import * as MVSUtil from '@molstar/mvs/util';
import { KinemageExtension } from '@molstar/kinemage-extension/behavior';
import * as interactivity from '@molstar/plugin-extension/interactivity';
import * as loaders from '@molstar/plugin-extension/loaders';
import { ViewerPluginUIViewModel, ViewerPluginViewModel } from '@molstar/viewer/view-models';

export const ExtensionMap = {
  // Mol* built-in extensions
  mvs: PluginSpec.Behavior(MolViewSpec),
  backgrounds: PluginSpec.Behavior(Backgrounds),
  'debug-helpers': PluginSpec.Behavior(DebugHelpers),
  'model-export': PluginSpec.Behavior(ModelExport),
  'mp4-export': PluginSpec.Behavior(Mp4Export),
  'geo-export': PluginSpec.Behavior(GeometryExport),
  'zenodo-import': PluginSpec.Behavior(ZenodoImport),
  'wwpdb-chemical-component-dictionary': PluginSpec.Behavior(wwPDBChemicalComponentDictionary),
  kinemage: PluginSpec.Behavior(KinemageExtension),

  // 3rd party extensions
  'pdbe-structure-quality-report': PluginSpec.Behavior(PDBeStructureQualityReport),
  'dnatco-ntcs': PluginSpec.Behavior(DnatcoNtCs),
  'assembly-symmetry': PluginSpec.Behavior(AssemblySymmetry),
  'rcsb-validation-report': PluginSpec.Behavior(RCSBValidationReport),
  'anvil-membrane-orientation': PluginSpec.Behavior(ANVILMembraneOrientation),
  g3d: PluginSpec.Behavior(G3DFormat), // TODO: consider removing this for Mol* 6.0
  'ma-quality-assessment': PluginSpec.Behavior(MAQualityAssessment),
  'sb-ncbr-partial-charges': PluginSpec.Behavior(SbNcbrPartialCharges),
  tunnels: PluginSpec.Behavior(SbNcbrTunnels),
};

export const PluginExtensions = {
  wwPDBStructConn: wwPDBStructConnExtensionFunctions,
  mvs: {
    MVSData,
    createBuilder: MVSData.createBuilder,
    loadMVS,
    loadMVSData,
    util: {
      ...MVSUtil,
    },
  },
  modelArchive: {
    qualityAssessment: {
      config: MAQualityAssessmentConfig,
    },
  },
  plugin: {
    interactivity,
    loaders,
    models: {
      PluginViewModel: ViewerPluginViewModel,
      PluginUIViewModel: ViewerPluginUIViewModel,
    },
  },
};
