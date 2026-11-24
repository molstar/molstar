/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved view of a section in groups: its tracks flat in spec order, each naming its group, the groups in
// `trackGroups`, a hidden group's tracks left out, an authored summary marked, and the one walker.

import type { Node } from '@molstar/core/repr/tree';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { tracksOf } from '../groups.js';
import { resolveView } from '../view.js';
import { groupsSpec } from './groups-spec.js';

const sectionOf = (spec: Node) => resolveView(resolveSpec(spec), '/views/v')!.sections[0];

describe('resolved view: track groups', () => {
  const spec = groupsSpec();
  const section = sectionOf(spec);

  it('the spec checks clean', () => expect(checkSeqSpec(spec)).toEqual([]));

  it('keeps the tracks flat in spec order, each naming its group; a hidden track left out', () => {
    expect(section.tracks.map((t) => [t.address, t.group ?? null])).toEqual([
      ['/views/v/p53/letters', null],
      ['/views/v/p53/features/domains', '/views/v/p53/features'],
      ['/views/v/p53/features/motifs', '/views/v/p53/features'],
      ['/views/v/p53/variants/clinvar', '/views/v/p53/variants'],
      ['/views/v/p53/variants/gnomad', '/views/v/p53/variants'],
      ['/views/v/p53/variants/density', '/views/v/p53/variants'],
      ['/views/v/p53/signals/a', '/views/v/p53/signals'],
      ['/views/v/p53/signals/b', '/views/v/p53/signals'],
    ]);
    // every layer placed on the section's system, its address under its track's
    expect(section.tracks.every((t) => t.layers.every((l) => l.placed && l.address?.startsWith(t.address)))).toBe(true);
  });

  it('says the groups in order: label, opening state, authored summary (marked on its track), count', () => {
    expect(section.trackGroups).toEqual([
      { address: '/views/v/p53/features', label: 'Features', collapsed: false, count: 2 },
      {
        address: '/views/v/p53/variants',
        label: 'Variants',
        collapsed: true,
        summary: '/views/v/p53/variants/density',
        count: 2,
      },
      { address: '/views/v/p53/signals', label: 'Signals', collapsed: false, count: 2 },
    ]);
    expect(section.tracks.filter((t) => t.summaryOf).map((t) => [t.address, t.summaryOf])).toEqual([
      ['/views/v/p53/variants/density', '/views/v/p53/variants'],
    ]);
  });

  it('`summary: none` is said; a hidden group leaves its tracks out and is no group', () => {
    const none = sectionOf(groupsSpec({ variants: { summary: 'none' } }));
    expect(none.trackGroups![1]).toEqual({
      address: '/views/v/p53/variants',
      label: 'Variants',
      collapsed: true,
      none: true,
      count: 3,
    });
    const hidden = sectionOf(groupsSpec({ signals: { hidden: true } }));
    expect(hidden.tracks.some((t) => t.group === '/views/v/p53/signals')).toBe(false);
    expect(hidden.trackGroups!.map((g) => g.address)).toEqual(['/views/v/p53/features', '/views/v/p53/variants']);
  });

  it('a section without groups has no trackGroups and its tracks no group (a 0.10 spec reads unchanged)', () => {
    const b = groupsSpec();
    const sec = b.children![1].children![0].children![0];
    sec.children = sec.children!.filter((c) => c.kind === 'track');
    const flat = sectionOf({ ...b, params: { ...b.params, version: '0.10.0' } });
    expect(flat.trackGroups).toBeUndefined();
    expect(flat.tracks.map((t) => [t.address, t.group])).toEqual([['/views/v/p53/letters', undefined]]);
  });

  it('the walker: tracks through groups, hidden included; a newer minor’s container read through (EC-244)', () => {
    const sec = spec.children![1].children![0].children![0];
    const withFuture: Node = {
      ...sec,
      children: [...sec.children!, { kind: 'bundle', id: 'x', children: [{ kind: 'track', id: 'later' }] }],
    };
    const walked = tracksOf(withFuture, '/s');
    expect(walked.map((t) => [t.address, t.group?.address ?? null])).toEqual([
      ['/s/letters', null],
      ['/s/features/domains', '/s/features'],
      ['/s/features/motifs', '/s/features'],
      ['/s/variants/clinvar', '/s/variants'],
      ['/s/variants/gnomad', '/s/variants'],
      ['/s/variants/density', '/s/variants'],
      ['/s/signals/a', '/s/signals'],
      ['/s/signals/b', '/s/signals'],
      ['/s/signals/c', '/s/signals'],
      // its container is no group: the track is the section's, as `group` is to an older reader
      ['/s/x/later', null],
    ]);
  });
});
