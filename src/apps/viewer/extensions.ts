/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ANVILMembraneOrientation } from '../../extensions/anvil/behavior';
import { AssemblySymmetry } from '../../extensions/assembly-symmetry';
import { Backgrounds } from '../../extensions/backgrounds';
import { DnatcoNtCs } from '../../extensions/dnatco';
import { G3DFormat } from '../../extensions/g3d/format';
import { GeometryExport } from '../../extensions/geo-export';
import { MAQualityAssessment, MAQualityAssessmentConfig } from '../../extensions/model-archive/quality-assessment/behavior';
import { ModelExport } from '../../extensions/model-export';
import { Mp4Export } from '../../extensions/mp4-export';
import { loadMVS } from '../../extensions/mvs';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { loadMVSData } from '../../extensions/mvs/components/formats';
import { PDBeStructureQualityReport } from '../../extensions/pdbe';
import { RCSBValidationReport } from '../../extensions/rcsb';
import { SbNcbrPartialCharges, SbNcbrTunnels } from '../../extensions/sb-ncbr';
import { wwPDBChemicalComponentDictionary } from '../../extensions/wwpdb/ccd/behavior';
import { wwPDBStructConnExtensionFunctions } from '../../extensions/wwpdb/struct-conn';
import { ZenodoImport } from '../../extensions/zenodo';
import { PluginSpec } from '../../mol-plugin/spec';
import { MVSData } from '../../extensions/mvs/mvs-data';
import * as MVSUtil from '../../extensions/mvs/util';

export const ExtensionMap = {
    // Mol* built-in extensions
    'mvs': PluginSpec.Behavior(MolViewSpec),
    'backgrounds': PluginSpec.Behavior(Backgrounds),
    'model-export': PluginSpec.Behavior(ModelExport),
    'mp4-export': PluginSpec.Behavior(Mp4Export),
    'geo-export': PluginSpec.Behavior(GeometryExport),
    'zenodo-import': PluginSpec.Behavior(ZenodoImport),
    'wwpdb-chemical-component-dictionary': PluginSpec.Behavior(wwPDBChemicalComponentDictionary),

    // 3rd party extensions
    'pdbe-structure-quality-report': PluginSpec.Behavior(PDBeStructureQualityReport),
    'dnatco-ntcs': PluginSpec.Behavior(DnatcoNtCs),
    'assembly-symmetry': PluginSpec.Behavior(AssemblySymmetry),
    'rcsb-validation-report': PluginSpec.Behavior(RCSBValidationReport),
    'anvil-membrane-orientation': PluginSpec.Behavior(ANVILMembraneOrientation),
    'g3d': PluginSpec.Behavior(G3DFormat), // TODO: consider removing this for Mol* 6.0
    'ma-quality-assessment': PluginSpec.Behavior(MAQualityAssessment),
    'sb-ncbr-partial-charges': PluginSpec.Behavior(SbNcbrPartialCharges),
    'tunnels': PluginSpec.Behavior(SbNcbrTunnels),
};

export const PluginExtensions = {
    wwPDBStructConn: wwPDBStructConnExtensionFunctions,
    mvs: {
        MVSData,
        createBuilder: MVSData.createBuilder,
        loadMVS,
        loadMVSData,
        util: {
            ...MVSUtil
        }
    },
    modelArchive: {
        qualityAssessment: {
            config: MAQualityAssessmentConfig
        }
    }
};
