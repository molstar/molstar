/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section of a source's complex (_spec/README.md "Overlay", Opened sources; "Resolved view", A section of a source's
// complex): resolved over the opened model once its source holds the complex, pending before, reflowed by an update
// (a unit dropped leaves its place, one gained takes one), the complex dropped said once; one rebuild a node version.

import type { Node } from '@molstar/core/repr/tree';
import { type NodeEntry, nodeKey } from '@molstar/seqx-core/data';
import { structureNodes } from '@molstar/seqx-core/io/structure';
import { checkSeqSpec, createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { Source, zincFinger } from '../../io/_test/fixtures.js';
import { resolveSpec } from '../../resolve/model.js';
import { openedModel } from '../../resolve/opened.js';
import { buildOverlay, EmptyOverlay, nodeAddress } from '../../resolve/overlay.js';
import { opensComplexes } from '../blocks.js';
import { resolveView } from '../view.js';

const Complex = nodeAddress(Source, 'complex');
const View = '/views/v';
const Section = `${View}/c`;

function spec(blocks?: (s: any) => void): Node {
  const b = createSeqSpec();
  b.data.source('zf', { format: 'mmcif', uri: 'zf.cif' });
  b.views.view('v', (v) =>
    v.section('c', { of: Complex }, (s) => {
      blocks?.(s);
      s.track('t').layer('letters', { type: 'letters' });
    }),
  );
  return b.build();
}

const All = structureNodes(zincFinger(), { source: Source }).nodes;

/** The nodes with the complex's units narrowed (`keep`), or the complex gone (`null`). */
function nodesWith(keep: ((id: string) => boolean) | null) {
  return All.flatMap((n) => {
    if (n.key !== 'complex') return [n];
    if (keep === null) return [];
    const children = (n.node.children ?? []).filter((c) => c.kind !== 'unit' || keep(c.id!));
    // the links name the units: kept only while every unit is
    return [{ ...n, node: { ...n.node, children: children.filter((c) => c.kind !== 'links' || keep('B-1_555')) } }];
  });
}

function reader(nodes: readonly { key: string; node: Node }[], version: number) {
  const map = new Map<string, NodeEntry>();
  for (const n of nodes) map.set(nodeKey(Source, n.key), { source: Source, key: n.key, node: n.node, version });
  return { nodes: map, nodeVersion: { value: version } };
}

/** The view resolved over the opened model of `nodes` at a node version. */
const viewOf = (
  model: ReturnType<typeof resolveSpec>,
  nodes: readonly { key: string; node: Node }[],
  version: number,
) => resolveView(openedModel(model, buildOverlay(model, reader(nodes, version))), View);
const withoutB = () => nodesWith((id) => id !== 'B-1_555');

const blocks = (v: ReturnType<typeof resolveView>) =>
  v!.sections[0].blocks?.map((b) => [b.label, b.systems.map((s) => s.slice(Complex.length + 1))]);
const codes = (v: ReturnType<typeof resolveView>) => v!.issues.map((i) => [i.code, i.address]);

describe('a section of a source’s complex', () => {
  it('the spec is valid, and a section of a source’s complex is what opens the model', () => {
    expect(checkSeqSpec(spec())).toEqual([]);
    const model = resolveSpec(spec());
    expect(opensComplexes(model)).toBe(true);
    const plain = createSeqSpec();
    plain.data.sequence('s', { alphabet: 'protein', letters: 'ACD' });
    plain.views.view('v', (v) => v.section('s', { of: '/data/s' }));
    expect(opensComplexes(resolveSpec(plain.build()))).toBe(false);
  });

  it('pending before its source holds a node: no block, no issue, an axis of no length', () => {
    const model = resolveSpec(spec());
    for (const m of [model, openedModel(model, EmptyOverlay)]) {
      const v = resolveView(m, View);
      expect(codes(v)).toEqual([]);
      expect(v!.sections[0].blocks).toBeUndefined();
      expect(v!.sections[0].axis.ranges).toBeUndefined();
    }
  });

  it('resolves once held: a block a unit (copies merged, labeled by id), its letters placed on the unit', () => {
    const model = resolveSpec(spec());
    const v = viewOf(model, All, 1)!;
    expect(codes(v)).toEqual([]);
    expect(blocks(v)).toEqual([
      ['A-1_555 ×2', ['A-1_555', 'B-1_555']],
      ['ZINC ION', ['C-1_555']],
    ]);
    expect(v.sections[0].blocks![0].layers[0][0]).toMatchObject({ placed: true });
    expect(v.sections[0].axis).toMatchObject({ system: Section, length: 6, ranges: [0, 6] });
  });

  it('reflows on an update: a unit dropped leaves its place, gained back takes it again', () => {
    const model = resolveSpec(spec());
    const dropped = viewOf(model, withoutB(), 2);
    expect(codes(dropped)).toEqual([]);
    expect(blocks(dropped)).toEqual([
      ['A-1_555', ['A-1_555']],
      ['ZINC ION', ['C-1_555']],
    ]);
    expect(blocks(viewOf(model, All, 3))).toEqual(blocks(viewOf(model, All, 1)));
  });

  it('an override of a unit an update drops goes with it: the rest stay derived', () => {
    const model = resolveSpec(spec((s) => s.block('B-1_555', { form: 'residues', label: 'B' })));
    const v = viewOf(model, withoutB(), 2);
    expect(codes(v)).toEqual([]);
    expect(blocks(v)).toEqual([
      ['A-1_555', ['A-1_555']],
      ['ZINC ION', ['C-1_555']],
    ]);
  });

  it('a complex its source no longer holds is said once, at the section; nothing is placed', () => {
    const model = resolveSpec(spec());
    const v = viewOf(model, nodesWith(null), 4);
    expect(codes(v)).toEqual([['seq.view.section.complex', Section]]);
    expect(v!.sections[0].blocks).toBeUndefined();
  });

  it('one opened model and one resolved view a node version: the same overlay gives the same objects', () => {
    const model = resolveSpec(spec());
    const r = reader(All, 1);
    const overlay = buildOverlay(model, r);
    const m = openedModel(model, overlay);
    expect(openedModel(model, overlay)).toBe(m);
    expect(resolveView(m, View)).toBe(resolveView(openedModel(model, overlay), View));
    const next = openedModel(model, buildOverlay(model, reader(All, 2)));
    expect(next).not.toBe(m);
    expect(next.hash).not.toBe(m.hash);
    // the static model's own reads are kept: its spec, its sources, its paths between its own systems
    expect([next.spec, next.sources, next.maps]).toEqual([model.spec, model.sources, model.maps]);
  });
});
