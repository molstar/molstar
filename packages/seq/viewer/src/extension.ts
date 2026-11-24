/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineExtension } from '@molstar/seqx-core/extension';
import type { SeqViewerExtension, SeqViewerSpec } from './viewer-spec.js';

// @molstar/seqx-viewer/extension (this file, through the package's `./*` export): a viewer extension's descriptor
// checked and completed where it is defined (_spec/README.md "Composition", Extensions), so a viewer composed with none
// bundles none of it.

/**
 * A viewer extension's descriptor: the core's checks and tables (`defineExtension`), and `composeViewer`, which adds its
 * viewer parts to a viewer's spec: its header entries to the registry, its behaviors after the spec's (its commands as a
 * behavior's, `<id>:commands`), its bindings a layer over the spec's, naming the viewer's actions or its own (`<id>:…`).
 */
export function defineViewerExtension<E extends SeqViewerExtension>(extension: E): E {
  const e = defineExtension(extension, true);
  const { id, commands, bindings = [], headerActions = [] } = e;
  return {
    ...e,
    composeViewer: (spec: SeqViewerSpec): SeqViewerSpec => ({
      ...spec,
      registry: [
        ...(spec.registry ?? []),
        ...headerActions.map((entry) => ({ kind: 'header-action' as const, entry })),
      ],
      behaviors: [
        ...spec.behaviors,
        ...(e.behaviors ?? []),
        ...(commands?.length ? [{ provider: { name: `${id}:commands`, create: () => ({ commands }) } }] : []),
      ],
      bindings: {
        ...spec.bindings,
        overrides: [...(spec.bindings?.overrides ?? []), bindings],
        actions: [
          ...(spec.bindings?.actions ?? []),
          ...bindings.flatMap((b) => (b.action?.startsWith(`${id}:`) ? [b.action] : [])),
        ],
      },
    }),
  };
}
