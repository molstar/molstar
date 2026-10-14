/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-a: dynamic nodes in the format (seqspec:format#nodes-and-refs, seqspec:sources#dynamic-nodes). A streamed
// source's derived refs are its dynamic nodes, valid in a document before they arrive (another source's are not); a
// data-mode resolver accepts them under the data's streamed sources; and `loadStreamNode` reads a node as a document
// node of its kind: checked (an unknown param, a ref that names nothing held, a bad key each reported and the node left
// out), its ref derived.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { type Issue, KeyKind, type KeyTarget, SeqSets } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  compile,
  indexDoc,
  type LoadResult,
  loadStreamNode,
  resolveView,
  type StageContext,
  type StreamNode,
  supportedAll,
  validate,
  withDynamic,
} from '../../src/format/index.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
const codes = (issues: readonly Issue[]) => issues.map((i) => i.code);

describe("a streamed source's derived refs are its dynamic nodes", () => {
  const doc = (stream: boolean, data: string) => {
    const b = createSeqViewSpec({ title: 't' });
    b.data.sequence({ ref: 'chr', value: 'ACGTACGTAC' });
    if (stream) b.data.source({ ref: 'src', stream: 's', data_kind: 'values', coordinate_system: 'chr' });
    else b.data.source({ ref: 'src', format: 'fasta', value: '>x\nACGT\n' });
    b.view({ ref: 'v', label: 'v' }).track().ranges({ data });
    return b.getDocument();
  };

  it("valid in a document before they arrive; another source's derived refs name nothing", () => {
    expect(codes(validate(doc(true, 'src/peaks-1'), reg))).toEqual([]);
    expect(codes(validate(doc(false, 'src/peaks-1'), reg))).toEqual(['svs.ref.wrong_kind']);
  });

  it("in data mode, under the data's streamed sources only", () => {
    const { set } = SeqSets.build({ molecules: [{ key: 'chr', residues: 'ACGT', ownUnit: true }] });
    const keys = withDynamic(set, new Set(['src']));
    const t: KeyTarget = { kind: -1, id: -1, row: -1 };
    expect([keys.dynamic?.('src/peaks-1'), keys.dynamic?.('other/x'), keys.dynamic?.('chr')]).toEqual([
      true,
      false,
      false,
    ]);
    expect(keys.byName('chr', t) && t.kind).toBe(KeyKind.molecule);
    const views = (data: string) => {
      const b = createSeqViewSpec({ title: 't' });
      b.view({ ref: 'v', label: 'v' }).track().ranges({ data });
      return b.getDocument();
    };
    expect(codes(validate(views('src/peaks-1'), reg, { set: keys }))).toEqual([]);
    expect(codes(validate(views('other/x'), reg, { set: keys }))).toEqual(['svs.ref.unresolved']);
    expect(withDynamic(set, new Set())).toBe(set);
  });
});

describe('a layer whose data waits for dynamic nodes', () => {
  async function issues(data: string | string[]) {
    const b = createSeqViewSpec({ title: 't' });
    const p = b.data.sequence({ ref: 'p', value: 'MKVLAAGG' });
    b.data.source({ ref: 'genes', stream: 'genes', data_kind: 'items' });
    b.data.source({ ref: 'peaks', stream: 'peaks', data_kind: 'items' });
    const v = b.view({ ref: 'v', label: 'v' });
    v.block({ sequence: p });
    v.track({ ref: 't' }).ranges({ data });
    const d = b.getDocument();
    const { set } = await compile(d, reg, { load: () => Promise.reject(new Error('no files')), stream: () => null });
    return resolveView(indexDoc(d, reg), set, 'v', reg, supportedAll(reg)).issues;
  }

  it('is noted once per source it waits for, with a count and the first refs (svs.ref.pending, info)', async () => {
    const said = await issues(['genes/a', 'genes/b', 'genes/c', 'genes/d', 'genes/a', 'peaks/x']);
    expect(said.map((i) => [i.severity, i.code, i.message])).toEqual([
      [
        'info',
        'svs.ref.pending',
        '"genes/a", "genes/b", "genes/c" and 1 more name nothing "genes" has sent yet: shown when they arrive',
      ],
      ['info', 'svs.ref.pending', '"peaks/x" names nothing "peaks" has sent yet: shown when it arrives'],
    ]);
    expect(said[0].path).toMatch(/\/params\/data$/);
  });
});

describe('loadStreamNode', () => {
  const load = async (node: StreamNode, known: string[] = []) => {
    const issues: Issue[] = [];
    const ctx: StageContext = {
      signal: new AbortController().signal,
      registries: reg,
      set: null,
      sources: new Map(),
      unavailable: new Set(),
      source: () => null,
      report: (i) => issues.push(i),
    };
    const part: LoadResult | null = await loadStreamNode('src', node, (r) => known.includes(r), ctx, fixtureLoaders());
    return { part, issues };
  };

  it('reads a node as a document node of its kind, under its derived ref', async () => {
    const seq = await load({ kind: 'sequence', key: 'p1', params: { value: 'MKVL', label: 'P1' } });
    expect(seq.issues).toEqual([]);
    expect(seq.part?.molecules).toMatchObject([{ key: 'src/p1', residues: 'MKVL', label: 'P1', ownUnit: true }]);
    const ann = await load(
      {
        kind: 'annotation',
        key: 'a/1',
        params: { geometry: 'ranges', target: 'src/p1/unit', value: { start: [0], stop: [2] } },
      },
      ['src/p1/unit'],
    );
    expect(ann.issues).toEqual([]);
    expect(ann.part?.annotations?.[0]).toMatchObject({ key: 'src/a/1', geometry: 'ranges', target: 'src/p1/unit' });
  });

  it('leaves a node out with its issue: an unknown param, a ref to nothing held, a bad key, an unknown kind', async () => {
    const typo = await load({ kind: 'sequence', key: 'p', params: { valeu: 'MKVL' } });
    expect([typo.part, codes(typo.issues)]).toEqual([null, ['svs.params.unknown', 'svs.params.exclusive']]);
    expect(typo.issues[0].message).toMatch(/^src: /);
    const lost = await load({ kind: 'unit', key: 'u', params: { sequence: 'src/p9' } });
    expect([lost.part, codes(lost.issues)]).toEqual([null, ['svs.stream.ref']]);
    expect(lost.issues[0]).toMatchObject({ severity: 'warning', path: '/nodes/u/params/sequence' });
    const key = await load({ kind: 'sequence', key: 'bad key', params: { value: 'MK' } });
    expect([key.part, codes(key.issues)]).toEqual([null, ['svs.ref.syntax']]);
    const kind = await load({ kind: 'layer' as never, key: 'l', params: {} });
    expect([kind.part, codes(kind.issues)]).toEqual([null, ['svs.node.kind']]);
    // A nested source streams (ST1b-b): one without a stream, or its data kind, is left out with its issue.
    const file = await load({ kind: 'source', key: 's', params: { format: 'fasta', value: '>x\nMK\n' } });
    expect([file.part, codes(file.issues)]).toEqual([null, ['svs.params.missing']]);
  });
});
