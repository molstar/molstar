/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the graph and document tests: a small document (inline sources, no network), stub renderers and
// themes that count their calls (declaring `content: false` as the kernel's will, seqspec:viewer#updates-and-status),
// and loaders whose requests a test resolves or fails by hand.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  KernelRegistry,
  type Loaders,
  type LoadRequest,
  type Patch,
  type SeqViewSpecDoc,
} from '@molstar/seq-core/format';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  type LayerRenderer,
  type SeqViewerContext,
  type SeqViewerPorts,
  type ThemeProvider,
} from '../src/index.js';

export const KERNEL_TYPES = ['letters', 'ruler', 'ranges', 'points', 'values', 'rows', 'color', 'pairs'];
/** The kernel renderers and themes that never read residues or alignment content. */
const NO_CONTENT = new Set(['ranges', 'points', 'variants', 'secondary', 'values', 'color', 'pairs', 'matrix']);

/** `DefaultRegistry` without its renderers and themes, for a test's own in their place. */
export const kernelFree = () => ({ ...DefaultRegistry, renderers: undefined, themes: undefined });

/** Renderers for the kernel types that draw linear and wrapped, counting prepares by layer ref. */
export function stubs() {
  const prepared = new Map<string, number>();
  const disposed: unknown[] = [];
  const renderers: LayerRenderer[] = KERNEL_TYPES.map((name) => ({
    name,
    geometries: ['linear', 'wrapped'],
    content: !NO_CONTENT.has(name),
    prepare: (input) => {
      prepared.set(input.layer.ref, (prepared.get(input.layer.ref) ?? 0) + 1);
      return { ref: input.layer.ref, rows: input.layer.placements.length };
    },
    dispose: (s) => void disposed.push(s),
  }));
  const created: string[] = [];
  // Every kernel theme, for the kinds its format entry serves; only `rules` precomputes from residues.
  const themes: ThemeProvider[] = (KernelRegistry.themes ?? []).map((t) => ({
    name: t.name,
    label: t.name,
    category: 'test',
    granularity: 0,
    kinds: Object.keys(t.params) as ThemeProvider['kinds'],
    content: t.name === 'rules',
    isApplicable: () => true,
    create: (tc, _params, channel) => {
      created.push(`${t.name}:${tc.layer.ref}`);
      return { kind: channel.kind, name: t.name, dark: tc.dark } as never;
    },
  }));
  return { prepared, disposed, created, entry: { renderers, themes } };
}

export function viewer(ports: Partial<SeqViewerPorts> = {}) {
  const s = stubs();
  const spec = DefaultSeqViewerSpec();
  const ctx = createSeqViewer({ ...spec, registry: [kernelFree(), s.entry] }, ports);
  return { ctx, ...s };
}

/** Two views over one sequence read from an inline FASTA, a numbering, an annotation, and a column statistic. */
export function smallDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const src = b.data.source({ ref: 'pasted', value: '>a\nMKVLATGHEK\n>b\nMKVLA\n', format: 'fasta' });
  b.data
    .sequence({ ref: 'a', source_ref: src, record: 'a' })
    .numbering({ ref: 'a.auth', name: 'auth', method: 'offset', first: 101 });
  b.data.sequence({ ref: 'b', source_ref: src, record: 'b' });
  b.data.annotation({
    ref: 'dom',
    geometry: 'ranges',
    target: 'a',
    value: { start: [0, 4], stop: [3, 8], category: ['x', 'y'] },
  });
  b.data.annotation({ ref: 'site', geometry: 'ranges', target: 'b', value: { start: [1], stop: [2] } });
  const v = b.view({ ref: 'v', label: 'V' });
  v.track({ ref: 'v.ruler' }).ruler({ ref: 'v.ruler.l' });
  v.track({ ref: 'v.letters' }).letters({ ref: 'v.letters.l' });
  v.track({ ref: 'v.dom' }).ranges({ ref: 'v.dom.l', data: 'dom' });
  v.track({ ref: 'v.site' }).ranges({ ref: 'v.site.l', data: 'site' });
  const w = b.view({ ref: 'w', label: 'W' });
  w.block({ sequence: 'b' });
  w.track({ ref: 'w.letters' }).letters({ ref: 'w.letters.l' });
  return b.getDocument();
}

/** Two inline sequences in one section: letters, and features on `y` (`content: false`); for the content keys. */
export function residuesDoc(x = 'MKVLATGHEK') {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'x', value: x }).numbering({ ref: 'x.auth', name: 'auth', method: 'offset', first: 101 });
  b.data.sequence({ ref: 'y', value: 'MKVLA' });
  b.data.annotation({ ref: 'feat', geometry: 'ranges', target: 'y', value: { start: [0], stop: [3] } });
  const v = b.view({ ref: 'v', label: 'V' });
  v.track({ ref: 'v.letters' }).letters({ ref: 'v.letters.l' });
  v.track({ ref: 'v.feat' }).ranges({ ref: 'v.feat.l', data: 'feat' });
  return b.getDocument();
}

/** `p` placed on `q`'s block through a segments mapping. */
export function mappedDoc() {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'p', value: 'MKVLAT' });
  b.data.sequence({ ref: 'q', value: 'MKVLATGH' });
  b.data.mapping({ ref: 'm', from: 'p', to: 'q', method: 'segments', segments: [[0, 0, 6]] } as never);
  b.data.annotation({ ref: 'site', geometry: 'ranges', target: 'p', value: { start: [0], stop: [2] } });
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ sequence: 'q' });
  v.track({ ref: 'v.site' }).ranges({ ref: 'v.site.l', data: 'site' });
  return b.getDocument();
}

/** `smallDoc` with `dom`'s second range ending at 9: one node's params changed. */
export function domMoved(): SeqViewSpecDoc {
  const d = smallDoc();
  const dom = d.root.children?.[0].children?.find((c) => c.ref === 'dom');
  if (dom) dom.params = { ...dom.params, value: { start: [0, 4], stop: [2, 9], category: ['x', 'y'] } };
  return d;
}

export const patch = (ops: Patch['ops'], extra: Partial<Patch> = {}): Patch => ({
  format: 'seq_view_spec_patch',
  ops,
  ...extra,
});

/** Loaders a test drives: each request waits until `resolve(uri, text)` or `fail(uri)`; they ignore aborts. */
export function manualLoaders() {
  const waiting = new Map<string, { resolve(b: Uint8Array): void; reject(e: Error): void }>();
  const requests: string[] = [];
  const loaders: Loaders = {
    load: (r: LoadRequest) =>
      new Promise((resolve, reject) => {
        const uri = String(r.uri);
        requests.push(uri);
        waiting.set(uri, { resolve, reject });
      }),
  };
  return {
    loaders,
    requests,
    resolve: (uri: string, text: string) => waiting.get(uri)?.resolve(new TextEncoder().encode(text)),
    fail: (uri: string) => waiting.get(uri)?.reject(new Error(`${uri}: 404`)),
  };
}

export const refsRan = (ctx: SeqViewerContext) => [...ctx.graph.ran].sort();
