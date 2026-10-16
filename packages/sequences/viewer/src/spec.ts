/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer spec (seqspec:viewer#the-viewer-spec), in the shape of Mol* v6's PluginSpec: registry entries, behaviors,
// config. The default composition is `default-spec.ts`, which nothing in the viewer imports: a host that lists its own
// parts bundles only those.

import type { SvsRegistryEntry } from '@molstar/seq-core/format';
import type { Severity } from '@molstar/seq-core/model';
import type { SeqBehaviorProvider } from './behaviors.js';
import type { SeqViewerConfigItem } from './config.js';
import type { RenderProps } from './layout.js';
import type { SeqViewerRegistryEntry } from './registry.js';
export interface SeqViewerSpec {
  /** Registered by createSeqViewer, in order, before any behavior is created. */
  readonly registry?: readonly SeqViewerRegistryEntry[];
  readonly behaviors: readonly SeqViewerSpec.Behavior[];
  readonly config?: readonly (readonly [SeqViewerConfigItem, unknown])[];
  /** The analog of Mol*'s canvas3d props: the initial `ctx.renderProps` over the defaults. */
  readonly render?: Partial<RenderProps>;
}
export declare namespace SeqViewerSpec {
  interface Behavior {
    readonly behavior: SeqBehaviorProvider<any>;
    readonly defaultParams?: object;
  }
}
export const SeqViewerSpec = {
  Behavior<P extends object>(behavior: SeqBehaviorProvider<P>, defaultParams?: Partial<P>): SeqViewerSpec.Behavior {
    return { behavior, defaultParams };
  },
};

/**
 * The brand of the viewer's own format entries, the only ones `register` lets name `seqv`. A registered symbol, not the
 * objects' identity: a page holding two copies of the viewer (an HMR update, a host that bundles it twice) makes its
 * spec with one and its context with the other. Not enumerable, so a spread copy (an extension's) is not the viewer's.
 */
const OWN = Symbol.for('seqstar.viewer.seqv');
/** Marks an entry as the viewer's own; pure (it marks the entry it is given), so entries no spec registers are left out. */
export const viewerOwn = /* @__NO_SIDE_EFFECTS__ */ (e: SvsRegistryEntry): SvsRegistryEntry =>
  Object.defineProperty(e, OWN, { value: true });
/** Whether `e` is one of the viewer's own format entries, from any copy of the viewer. */
export const isViewerOwn = (e: SvsRegistryEntry): boolean => (e as { [OWN]?: boolean })[OWN] === true;

/**
 * The viewer's own issue codes (seqspec:viewer#registration), declared in the format's half so validation knows them.
 */
const codes = (severity: Severity, names: string) => names.split(' ').map((c) => ({ code: `seqv.${c}`, severity }));
export const ViewerIssueCodes: SvsRegistryEntry = viewerOwn({
  issueCodes: [
    ...codes('warning', 'patch.dropped patch.skipped history.refused port.refused issues.truncated ui.slot'),
    // Ring labels that do not fit are hidden: that depends on the width, not the input (seqspec:format#validation).
    ...codes('info', 'render.labels'),
    ...codes('error', 'renderer.error behavior.error command.error state.error target.range'),
  ],
});

/** The viewer's own layer type (seqspec:viewer#extension-points): other tracks' overviews over the whole axis. */
