/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseJson, printJson } from '@molstar/core/repr/codecs/json';
import { validateAndNormalize } from '@molstar/core/repr/schema';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { createSeqSpec, type GroupBuilder } from '../builder.js';
import { checkSeqSpec } from '../check.js';
import { FlattenedKinds } from '../group-language.js';
import { SeqLanguage, SeqLanguageVersion } from '../language.js';

/** A section of three tracks, two in a group, built through the builder. */
function grouped(group: Record<string, unknown> = {}) {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
  let g!: GroupBuilder;
  b.views.view('v', (v) => {
    const sec = v.section('sec', { of: s });
    sec.track('letters', {}, (t) => t.layer('l', { type: 'letters' }));
    g = sec.group('g', { label: 'Group', ...group });
    g.track('a', {}, (t) => t.layer('l', { type: 'letters' }));
    g.track('b', {}, (t) => t.layer('l', { type: 'ruler' }));
  });
  return { b, g };
}

const codes = (n: Node) => checkSeqSpec(n).map((i) => [i.code, i.severity]);

describe('the group kind (language 0.11.0)', () => {
  it('is the language 0.11.0’s: a section’s child, tracks its children; the builder writes it', () => {
    expect(SeqLanguageVersion).toBe('0.11.0');
    expect(SeqLanguage.kinds.group.parents).toEqual(['section']);
    expect(SeqLanguage.kinds.track.parents).toEqual(['section', 'group']);
    const { b, g } = grouped({ collapsed: true, summary: 'b' });
    expect(g.address).toBe('/views/v/sec/g');
    expect(String(g)).toBe('/views/v/sec/g');
    const spec = b.build();
    const sec = spec.children![1].children![0].children![0];
    expect(sec.children!.map((c) => [c.kind, c.id])).toEqual([
      ['track', 'letters'],
      ['group', 'g'],
    ]);
    expect(sec.children![1].children!.map((c) => c.id)).toEqual(['a', 'b']);
    expect(sec.children![1].params).toMatchObject({ label: 'Group', collapsed: true, summary: 'b' });
    expect(checkSeqSpec(spec)).toEqual([]);
  });

  it('`summary` is `auto`, `none` or one of its own tracks (`seq.lang.group.summary`)', () => {
    for (const summary of ['auto', 'none', 'a', 'b']) expect(codes(grouped({ summary }).b.build())).toEqual([]);
    expect(codes(grouped({ summary: 'letters' }).b.build())).toEqual([['seq.lang.group.summary', 'error']]);
    expect(codes(grouped({ summary: 'nothing' }).b.build())).toEqual([['seq.lang.group.summary', 'error']]);
  });

  it('is one level: a group in a group is refused as any parent the description refuses', () => {
    const spec = grouped().b.build();
    const g = spec.children![1].children![0].children![0].children![1];
    g.children!.push({ kind: 'group', id: 'inner', children: [] });
    expect(codes(spec)).toEqual([['repr.schema.parent.invalid', 'error']]);
  });

  it('a 0.10 spec (no groups) reads unchanged: valid, and its JSON round trip the same tree', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.view('v', (v) => v.section('sec', { of: s }).track('t', {}, (t) => t.layer('l', { type: 'letters' })));
    const at10: Node = { ...b.build(), params: { version: '0.10.0' } };
    expect(checkSeqSpec(at10)).toEqual([]);
    const read = parseJson(printJson(at10, { language: SeqLanguage }));
    expect(read.issues).toEqual([]);
    const back = validateAndNormalize(SeqLanguage, read.node as Node);
    expect(back.issues).toEqual([]);
    expect(hashNode(back.node as Node)).toBe(hashNode(at10));
  });

  it('an older reader flattens `group` (EC-244); this reader keeps a newer minor’s container with a warning', () => {
    expect(FlattenedKinds).toEqual({ group: '0.11.0' });
    // a spec of 0.10 holding a group: its writer told that a 0.10 reader flattens it (once a kind)
    const older = grouped().b.build();
    expect(codes({ ...older, params: { ...older.params, version: '0.10.0' } })).toEqual([
      ['seq.lang.version.kind', 'warning'],
    ]);
    expect(codes(grouped().b.build())).toEqual([]);
    // a 0.12 container this reader lacks, holding a track: one warning (its children are read through it, Resolved view)
    const spec = grouped().b.build();
    const sec = spec.children![1].children![0].children![0];
    sec.children!.push({ kind: 'bundle', id: 'x', children: [{ kind: 'track', id: 't', children: [] }] });
    expect(codes({ ...spec, params: { version: '0.12.0' } })).toEqual([['repr.schema.kind.unknown', 'warning']]);
    // at the reader's own minor it is an error
    expect(codes(spec)).toEqual([['repr.schema.kind.unknown', 'error']]);
  });
});
