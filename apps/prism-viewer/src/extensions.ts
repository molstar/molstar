/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AnvilExtension } from '@molstar/anvil-extension/prism-viewer';
import { AssemblySymmetryExtension } from '@molstar/assembly-symmetry-extension/prism-viewer';
import { CommandLanguageExtension } from '@molstar/command-language-extension/prism-viewer';
import { GeoExportExtension } from '@molstar/geo-export-extension/prism-viewer';
import { KinemageExtension } from '@molstar/kinemage-extension/prism-viewer';
import { ModelArchiveExtension } from '@molstar/model-archive-extension/prism-viewer';
import { ModelExportExtension } from '@molstar/model-export-extension/prism-viewer';
import type { LazyViewerExtension, ViewerExtension } from '@molstar/prism-viewer/spec';
import { LazyVideoExport } from '@molstar/video-export-extension/prism-viewer-stub';
import { LazyZenodo } from '@molstar/zenodo-extension/prism-viewer-stub';
import { DataOnlyExtensions } from './data-extensions.js';

/**
 * Mol*'s extensions the app ships, after the viewer's built-in features (viewer-spec.md §5, app.md §1):
 * those with operations and views, the heavy ones lazy (their stubs: the chunk loads on first use), then the data-only
 * ones (their plugin parts and annotations).
 */
export const MolstarViewerExtensions: readonly (ViewerExtension | LazyViewerExtension)[] = [
  ModelArchiveExtension,
  ModelExportExtension,
  GeoExportExtension,
  KinemageExtension,
  AssemblySymmetryExtension,
  AnvilExtension,
  CommandLanguageExtension,
  LazyVideoExport,
  LazyZenodo,
  ...DataOnlyExtensions,
];
