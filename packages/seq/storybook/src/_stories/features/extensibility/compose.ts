/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue, Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import type { SeqViewerSpec } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DemoExtension } from './demo-extension.js';

// The Extensibility story's composition (shared by the story and its headless test): one spec naming the example
// extension's layer type (with its fallback), its key on a built-in section and its scheme on a built-in layer; a
// viewer composed without the extension, which `extend` registers at run time and its undo takes away.

/** The first 60 residues of p53. */
export const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP';

/** The sites' starts (0-based): S6, S15, S33, S46 of p53's N terminus, each a serine of `Letters`. */
export const SiteStarts = [5, 14, 32, 45] as const;

export const SitesLayer = '/views/main/protein/sites/lollipops';
export const DomainsLayer = '/views/main/protein/domains/domains';

/** The spec: the lollipops over four sites (`ranges` where a reader lacks them), the domains in the warm scheme. */
export function extensibilitySpec(): Node {
  const b = createSeqSpec({ title: 'An extension, toggled', extensions: [DemoExtension] });
  const p = b.data.sequence('p53', { alphabet: 'protein', letters: Letters });
  const sites = b.data.ranges('sites', {
    target: p,
    rows: {
      start: new Int32Array(SiteStarts),
      end: new Int32Array(SiteStarts.map((x) => x + 1)),
      id: ['S6', 'S15', 'S33', 'S46'],
      label: ['S6', 'S15', 'S33', 'S46'],
    },
  });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: {
      start: new Int32Array([0, 39]),
      end: new Int32Array([39, 60]),
      id: ['TAD1', 'TAD2'],
      label: ['TAD1', 'TAD2'],
      category: ['transactivation', 'proline-rich'],
    },
  });
  b.views.view('main', (v) => {
    // the extension's key on a built-in section: read by its lollipops, carried by a reader without it
    const s = v.section('protein', { of: p, 'demo:stems': 'short' } as never);
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('sites', { label: 'Sites' }).layer('lollipops', {
      type: 'demo:lollipops',
      fallback: 'ranges',
      data: sites,
      head: 5,
    } as never);
    s.track('domains', { label: 'Domains' }).layer('domains', {
      type: 'ranges',
      data: domains,
      color: theme.scheme('demo:warm', { field: 'category' }),
    });
  });
  return b.build();
}

/** The default viewer's spec, headless when `ui` is undefined: the extension is not listed (`extend` adds it). */
export function extensibilityComposition(ui: boolean): SeqViewerSpec {
  const base = DefaultSeqViewerSpec();
  return ui ? base : { ...base, ui: undefined };
}

/** What the view and its frames say, a line each: the fallback and the scheme it lacks while the extension is off. */
export function issueLines(viewer: SeqViewer): string[] {
  const all: Issue[] = [...(viewer.view.resolved.value?.issues ?? []), ...viewer.layoutIssues.value];
  return all.length
    ? all.map((i) => `${i.severity} ${i.code}${i.address ? ` at ${i.address}` : ''}: ${i.message}`)
    : ['No issues: every layer is drawn as its own type.'];
}

export { DemoExtension };
