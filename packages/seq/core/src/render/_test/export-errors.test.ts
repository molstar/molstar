/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A spec with an error exported (_spec/README.md, Render, `renderSvg`): drawn as the canvas draws it, its errors in the
// picture's `<desc>` (escaped, the first ten and a count); a clean spec unchanged; a spec past the bounds blank.

import type { Issue, Node } from '@molstar/core/repr/tree';
import { createSeqSpec, DefaultSpecBounds, SeqLanguageVersion } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { withIssues } from '../export.js';
import { DrawLists, drawListsToSvg } from '../index.js';
import { renderSvg } from '../render-svg.js';

const size = { width: 400, height: 60 } as const;

function spec(data?: string): Node {
  const b = createSeqSpec({ title: 'errors' });
  const seq = b.data.sequence('hello', { alphabet: 'protein', letters: 'MKTAYIAKQRQISFVKSHFSRQ' });
  b.views.view('main', (v) => {
    const s = v.section('seq', { of: seq });
    s.track('letters').layer('letters', { type: 'letters' });
    if (data) s.track('hits').layer('hits', { type: 'ranges', data: data as never });
  });
  return b.build();
}

const error = (message: string, n = 0): Issue => ({
  code: `seq.test.e${n}`,
  severity: 'error',
  message,
  address: `/x/${n}`,
});
const accessible = { title: 't', description: 'A view', tracks: [], highlights: [] };
const descOf = (svg: string) => /<desc[^>]*>([^<]*)<\/desc>/.exec(svg)?.[1] ?? '';

describe('exportSvg of a spec with an error', () => {
  it('draws what resolves and says the error in the description; a clean spec has no such text', async () => {
    const bad = await renderSvg(spec('/data/helo'), { size });
    expect(bad.lists.height).toBeGreaterThan(0);
    expect(descOf(bad.svg)).toContain('The spec has an error');
    const clean = await renderSvg(spec(), { size });
    expect(clean.svg).not.toContain('The spec has');
    expect(withIssues(accessible, [])).toBe(accessible);
  });

  it('a newer minor’s unknown kind is a warning: drawn as what is known, no error in the description (EC-235)', async () => {
    const future = spec();
    const newer: Node = {
      ...future,
      params: { ...future.params, version: '0.11.0' },
      children: [...future.children!, { kind: 'from-the-future', id: 'x' }],
    };
    const r = await renderSvg(newer, { size });
    expect(r.issues.map((i) => [i.code, i.severity])).toEqual([['repr.schema.kind.unknown', 'warning']]);
    expect(r.lists.height).toBeGreaterThan(0);
    expect(descOf(r.svg)).not.toContain('The spec has an error');
    // at the reader's own minor the same kind is an error
    const own = await renderSvg({ ...newer, params: { ...newer.params, version: SeqLanguageVersion } }, { size });
    expect(own.issues.map((i) => i.severity)).toEqual(['error']);
  });

  it('escapes what an issue says (`<`, `&`, `]]>`) in the description', () => {
    const lists = new DrawLists();
    lists.width = 10;
    lists.height = 10;
    const accessibleWith = withIssues(accessible, [error('a <b> & ]]> c')]);
    const svg = drawListsToSvg(lists, { accessible: accessibleWith });
    expect(descOf(svg)).toContain('a &lt;b&gt; &amp; ]]&gt; c');
    expect(svg).not.toContain('<b>');
  });

  it('names the first ten errors and counts the rest', () => {
    const errors = Array.from({ length: 13 }, (_, i) => error(`m${i}`, i));
    const d = withIssues(accessible, errors).description;
    expect(d).toContain('seq.test.e9');
    expect(d).not.toContain('seq.test.e10');
    expect(d).toMatch(/ and 3 more$/);
  });

  it('a spec past the bounds is still a blank document with the bound’s issue', async () => {
    let node: Node = { kind: 'ext:note', id: 'c' };
    for (let i = 1; i < 5; i++) node = { kind: 'ext:note', id: 'c', children: [node] };
    const deep: Node = { kind: 'spec', params: { version: SeqLanguageVersion }, children: [node] };
    const r = await renderSvg(deep, { size, bounds: { ...DefaultSpecBounds, maxDepth: 3 } });
    expect(r.issues.map((i) => i.code)).toContain('seq.lang.bound.depth');
    expect(r.lists.passes.data.rects.n).toBe(0);
    expect(r.svg).not.toContain('<desc');
  });
});
