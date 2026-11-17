/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { type NodeEntry, nodeKey, SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { buildOverlay, EmptyOverlay, nodeAddress, resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { expandTemplate, templateNeeds, templateNodes } from '../templates.js';
import { resolveViews } from '../view.js';

// Track templates: the resolved view holds the template (no instance, whatever the store holds); expanded over an
// overlay, a track per selected node in its order, filtered, its data-less layers placed on the node.

const S = '/data/genes';

function setup(forEach: Record<string, unknown>, kind = 'ranges') {
  const b = createSeqSpec();
  const genes = b.data.source('genes', { format: 'seq', uri: 'genes.seq' });
  b.data.sequence('gene', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  b.data.sequence('gene2', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  b.views.view('main', (v) =>
    v
      .section('gene', { of: '/data/gene' })
      .track('t', { label: 'T', for_each: { source: genes, node_kind: kind, ...forEach } as never })
      .layer('marks', { type: 'ranges' }),
  );
  const model = resolveSpec(b.build(), { runtimes: [SeqRuntime] });
  return { model, view: resolveViews(model)[0] };
}

const ranges = (start: number, coding: boolean, target = '/data/gene'): Node => ({
  kind: 'ranges',
  params: { target, ...(coding ? { map: 'codon' } : {}), rows: { start: [start], end: [start + 9] } },
});
const source = (uri: string): Node => ({ kind: 'source', params: { format: 'seq', uri } });

function reader(nodes: [string, Node][]) {
  const map = new Map<string, NodeEntry>();
  for (const [key, node] of nodes) map.set(nodeKey(S, key), { source: S, key, node, version: 1 });
  return { nodes: map, nodeVersion: { value: 1 } };
}

const nodes = reader([
  ['b', ranges(50, true, '/data/gene2')],
  ['a', ranges(90, false)],
  ['c', ranges(10, true)],
  ['s1', source('z.seq')],
  ['s2', source('a.seq')],
  ['s3', source('m.seq')],
]);

describe('track templates', () => {
  it('stay a template in the resolved view, and list the depth of nesting they need', () => {
    const { view } = setup({ depth: 1 });
    const track = view.sections[0].tracks[0];
    expect(track.template).toEqual({ source: S, kind: 'ranges', depth: 1, order: 'key', filter: {} });
    expect(track.layers[0].placed).toBe(false);
    expect(view.issues).toEqual([]);
    expect(templateNeeds(view)).toEqual(new Map([[S, 1]]));
  });

  it('select nodes of their kind in key order, or by a param, descending with a minus', () => {
    const byKey = setup({});
    const o = buildOverlay(byKey.model, nodes);
    expect(templateNodes(byKey.view.sections[0].tracks[0], o).map((n) => n.key)).toEqual(['a', 'b', 'c']);
    const desc = setup({ order: '-key' });
    expect(templateNodes(desc.view.sections[0].tracks[0], o).map((n) => n.key)).toEqual(['c', 'b', 'a']);
    // by a param, ties by address
    const byTarget = setup({ order: 'target' });
    expect(templateNodes(byTarget.view.sections[0].tracks[0], o).map((n) => n.key)).toEqual(['a', 'c', 'b']);
    const byTargetDesc = setup({ order: '-target' });
    expect(templateNodes(byTargetDesc.view.sections[0].tracks[0], o).map((n) => n.key)).toEqual(['b', 'a', 'c']);
  });

  it('keep the nodes whose params match the filter', () => {
    const { model, view } = setup({ filter: { map: 'codon' }, order: '-key' });
    const o = buildOverlay(model, nodes);
    expect(templateNodes(view.sections[0].tracks[0], o).map((n) => n.key)).toEqual(['c', 'b']);
    const none = setup({ filter: { map: 'other' } });
    expect(templateNodes(none.view.sections[0].tracks[0], o)).toEqual([]);
  });

  it('expand into instance tracks whose layers draw their node, kept per overlay', () => {
    const { model, view } = setup({});
    const section = view.sections[0];
    const o = buildOverlay(model, nodes);
    const tracks = expandTemplate(model, o, section, section.tracks[0]);
    expect(tracks.map((t) => [t.label, t.template])).toEqual([
      ['T a', undefined],
      ['T b', undefined],
      ['T c', undefined],
    ]);
    const layer = tracks[0].layers[0];
    expect(tracks[0].address).toBe(`/views/main/gene/t/(node :at ${nodeAddress(S, 'a')})`);
    expect(layer).toMatchObject({
      address: `${tracks[0].address}/marks`,
      placed: true,
      data: { address: nodeAddress(S, 'a'), system: '/data/gene', shape: 'ranges' },
    });
    expect(expandTemplate(model, o, section, section.tracks[0])).toBe(tracks);
    expect(expandTemplate(model, EmptyOverlay, section, section.tracks[0])).toEqual([]);
    // at most `max` instances, the first in order
    expect(expandTemplate(model, o, section, section.tracks[0], 2).map((t) => t.label)).toEqual(['T a', 'T b']);
  });
});
