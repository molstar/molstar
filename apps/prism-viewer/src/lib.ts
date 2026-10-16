/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as LinearAlgebra3D from '@molstar/core/math/linear-algebra/3d';
import { DataLoci, EveryLoci, Loci } from '@molstar/model/model/loci';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import * as Structure from '@molstar/model/model/structure';
import { Volume } from '@molstar/model/model/volume';
import { MAQualityAssessmentConfig } from '@molstar/model-archive-extension/quality-assessment/behavior';
import { loadMVS } from '@molstar/mvs';
import { loadMVSData } from '@molstar/mvs/components/formats';
import * as MVSUtil from '@molstar/mvs/util';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { PluginConfig } from '@molstar/plugin/config';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultRegistry } from '@molstar/plugin/default-registry';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import { StateActions } from '@molstar/plugin/state/actions';
import { PluginStateObject, PluginStateTransform } from '@molstar/plugin/state/objects';
import * as interactivity from '@molstar/plugin-extension/interactivity';
import * as loaders from '@molstar/plugin-extension/loaders';
import { node } from '@molstar/prism-viewer/api';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import { MolScriptBuilder } from '@molstar/query-language/language/builder';
import { wwPDBStructConnExtensionFunctions } from '@molstar/wwpdb-extension/struct-conn';
import { FullPrismViewerSpec } from './full-spec.js';

/**
 * `molstar.lib` (app.md §1): Mol*'s library for scripts on the page, without a build step. The keys mirror v5's Viewer
 * (`apps/viewer/src/lib.ts`) where the Prism app has them; left out are the classic UI's (`PluginUIContext`,
 * `DefaultPluginUISpec`, the extensions' view models) and `StateTransforms` (the classic Viewer's 5.x carve-out,
 * `.v6/designs/architecture.md` §9.3; import a transformer from its module instead). `prism` holds the Prism viewer's
 * own entry points. The app loads it beside the viewer, in a chunk of its own; `@molstar/prism-viewer` does not
 * export it.
 */
export const lib = {
  structure: {
    ...Structure,
  },
  volume: {
    Volume,
  },
  shape: {
    Shape,
    ShapeGroup,
  },
  loci: {
    Loci,
    DataLoci,
    EveryLoci,
  },
  math: {
    LinearAlgebra: {
      ...LinearAlgebra3D,
    },
  },
  plugin: {
    PluginContext,
    PluginConfig,
    PluginBehavior,
    PluginSpec,
    PluginStateObject,
    PluginStateTransform,
    StateActions,
    DefaultPluginSpec,
    DefaultRegistry,
  },
  extensions: {
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
    },
  },
  prism: {
    createPrismViewer,
    DefaultPrismViewerSpec,
    FullPrismViewerSpec,
    node,
    MS: MolScriptBuilder,
  },
};

export type PrismViewerLib = typeof lib;
