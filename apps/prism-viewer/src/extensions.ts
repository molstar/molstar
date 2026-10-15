/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ModelArchiveExtension } from '@molstar/model-archive-extension/prism-viewer';
import type { ViewerExtension, LazyViewerExtension } from '@molstar/prism-viewer/spec';

/** Mol*'s extensions the app ships, after the viewer's built-in features (viewer-spec.md §5, plan phase 5). */
export const MolstarViewerExtensions: readonly (ViewerExtension | LazyViewerExtension)[] = [ModelArchiveExtension];
