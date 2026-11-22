/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The spec's tooltips headless (_spec/README.md "Render", Tooltips and Tooltip tables): a `ranges` table an annotation
// at its tooltip's address, read through the store and wanted beside the sections it decorates (and not otherwise);
// `tooltipsAt` by data and by range against a brute-force overlap; the cap and the count past it (EC-259); a table not
// yet come (EC-257); values escaped, the parse cleaned, Seq*'s own links alone (EC-256); MVS's schemas as `elements`
// forms, a row's other columns ignored (EC-252); addresses located by the caller.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import { type Location, rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { openModel } from '../export.js';
import { rowAddress, type TooltipOut, type TooltipPlace, tooltipsAt } from '../tooltips.js';
import { settle } from './helpers.js';

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP';
const P = '/data/p';
const Size = { size: { width: 900, height: 120 } };

type Build = ReturnType<typeof createSeqSpec>;

/** A protein with a domain layer; `tips` adds the tooltips. */
function spec(tips: (b: Build, p: string, domains: string) => void) {
  const b = createSeqSpec({ title: 'tooltips' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: { start: [2, 30], end: [12, 40], id: ['a', 'b'], label: ['A', 'B'], note: ['first **domain**', 'B*07:02'] },
  });
  b.views.view('v', (v) =>
    v.section('s', { of: p }, (s) => {
      s.track('letters').layer('l', { type: 'letters' });
      s.track('domains').layer('d', { type: 'ranges', data: domains });
    }),
  );
  tips(b, p, domains);
  return b.build();
}

const at = (position: number, extra: Partial<TooltipPlace> = {}): TooltipPlace => ({
  parts: [{ system: P, start: position, end: position + 1 }],
  ...extra,
});
const fresh = (): TooltipOut => ({ rows: [], more: 0, failed: [] });

describe('tooltip tables', () => {
  it('a `ranges` table is an annotation at its tooltip’s address, inline through the store', async () => {
    const s = await settle(
      spec((b, p) => b.views.tooltip('t', { target: p, rows: { start: [0], end: [5], tooltip: ['start'] } })),
      [SeqRuntime],
    );
    const a = s.model.annotations.get('/views/t');
    expect(a?.shape).toBe('ranges');
    expect(a?.access?.kind).toBe('inline');
    expect(s.view.tooltips?.map((t) => [t.address, t.schema])).toEqual([['/views/t', 'ranges']]);
  });

  it('its tiles are wanted beside the windows it decorates, and move with them; a view without one wants none', async () => {
    const tip = spec((b, p) => b.views.tooltip('t', { target: p, rows: { start: [0], end: [5], tooltip: ['x'] } }));
    const s = await settle(tip, [SeqRuntime]);
    const f = await s.frame({ ...Size, viewport: { start: 0, end: 20 } });
    const mine = f.wants.filter((w) => w.layer === '/views/t');
    expect(mine.length).toBeGreaterThan(0);
    expect(mine.every((w) => w.start <= 20 && w.end > 0)).toBe(true);
    const g = await s.frame({ ...Size, viewport: { start: 30, end: 50 } });
    expect(g.wants.filter((w) => w.layer === '/views/t').every((w) => w.end > 30)).toBe(true);
    const plain = await settle(
      spec(() => {}),
      [SeqRuntime],
    );
    const h = await plain.frame(Size);
    expect(h.wants.some((w) => w.layer?.startsWith('/views/') && !w.layer.includes('/v/'))).toBe(false);
  });
});

describe('tooltipsAt', () => {
  it('by range: every row covering the place, a row each, in position order; none elsewhere', async () => {
    const s = await settle(
      spec((b, p) =>
        b.views.tooltip('t', {
          target: p,
          rows: { start: [0, 3, 40], end: [10, 5, 41], tooltip: ['wide', 'narrow', 'far'] },
          label: 'Note',
        }),
      ),
      [SeqRuntime],
    );
    await s.frame(Size);
    const ctx = { model: s.model, view: s.view, data: s.store };
    const out = fresh();
    tooltipsAt(ctx, at(4), out);
    expect(out.rows).toEqual([
      { label: 'Note', value: 'wide', item: true },
      { label: 'Note', value: 'narrow', item: true },
    ]);
    const one = fresh();
    tooltipsAt(ctx, at(40), one);
    expect(one.rows).toEqual([{ label: 'Note', value: 'far' }]);
    const none = fresh();
    tooltipsAt(ctx, at(20), none);
    expect(none.rows).toEqual([]);
  });

  it('matches a brute-force overlap over random tables (a property)', async () => {
    let seed = 7;
    const rand = () => {
      seed = (seed * 16807) % 2147483647;
      return seed / 2147483647;
    };
    for (let round = 0; round < 3; round++) {
      const n = 60;
      const start = Array.from({ length: n }, () => Math.floor(rand() * Letters.length));
      const end = start.map((x) => Math.min(Letters.length, x + Math.floor(rand() * 8)));
      const tooltip = start.map((_, i) => `r${i}`);
      const s = await settle(
        spec((b, p) => b.views.tooltip('t', { target: p, rows: { start, end, tooltip }, format: 'text' })),
        [SeqRuntime],
      );
      await s.frame(Size);
      for (let q = 0; q < Letters.length; q += 3) {
        const out = fresh();
        tooltipsAt({ model: s.model, view: s.view, data: s.store }, at(q), out, 1000);
        // a row meets a position when it covers it (a zero-length row at its start)
        const brute = tooltip.filter((_, i) => start[i] <= q && Math.max(end[i], start[i] + 1) > q);
        expect(out.rows.map((r) => r.label).sort()).toEqual(brute.sort());
      }
    }
  });

  it('past the cap: the rows made and the rest counted (EC-259)', async () => {
    const n = 50;
    const s = await settle(
      spec((b, p) =>
        b.views.tooltip('t', {
          target: p,
          rows: {
            start: Array(n).fill(0),
            end: Array(n).fill(20),
            tooltip: Array.from({ length: n }, (_, i) => `n${i}`),
          },
        }),
      ),
      [SeqRuntime],
    );
    await s.frame(Size);
    const out = fresh();
    tooltipsAt({ model: s.model, view: s.view, data: s.store }, at(5), out, 3);
    expect(out.rows.map((r) => r.label)).toEqual(['n0', 'n1', 'n2']);
    expect(out.rows.every((r) => r.item)).toBe(true);
    expect(out.more).toBe(n - 3);
  });

  it('a table not yet come gives a loading row (EC-257)', async () => {
    const tip = spec((b, p) => b.views.tooltip('t', { target: p, rows: { start: [0], end: [5], tooltip: ['x'] } }));
    const store = new DataStore();
    const { model, view } = await openModel(tip, store, { runtimes: [SeqRuntime] });
    const out = fresh();
    tooltipsAt({ model, view: view!, data: store }, at(1), out);
    expect(out.rows).toEqual([{ label: 'loading', value: '', loading: true }]);
  });

  it('by data: the hit’s item’s own row, its fields read by name; another layer’s item none', async () => {
    const s = await settle(
      spec((b, _, domains) => b.views.tooltip('t', { on: domains, text: '{label}: {note}', label: 'Domain' })),
      [SeqRuntime],
    );
    await s.frame(Size);
    const layer = s.view.sections[0].tracks[1].layers[0];
    const letters = s.view.sections[0].tracks[0].layers[0];
    const ctx = { model: s.model, view: s.view, data: s.store };
    const out = fresh();
    tooltipsAt(ctx, at(4, { layer, item: { key: 'a', start: 2 } }), out);
    expect(out.rows).toHaveLength(1);
    expect(out.rows[0].label).toBe('Domain');
    expect(out.rows[0].value).toBe('A: first **domain**');
    const other = fresh();
    tooltipsAt(ctx, at(4, { layer: letters }), other);
    expect(other.rows).toEqual([]);
  });

  it('values escaped, the template’s Markdown read, cleaned after the parse, links Seq*’s alone (EC-256)', async () => {
    const s = await settle(
      spec((b, _, domains) => {
        b.views.tooltip('t', {
          on: domains,
          text: '**{note}** [see](!seq-focus=/refs/x) [away](https://x.org) &#x202E;',
        });
      }),
      [SeqRuntime],
    );
    await s.frame(Size);
    const layer = s.view.sections[0].tracks[1].layers[0];
    const out = fresh();
    tooltipsAt({ model: s.model, view: s.view, data: s.store }, at(31, { layer, item: { key: 'b', start: 30 } }), out);
    const [row] = out.rows;
    // the value's `*` stays text; the external link is its text; the bidi override decoded by the parse is removed
    expect(row.label).toBe('B*07:02 see away ');
    const rich = JSON.stringify(row.rich);
    expect(rich).toContain('"kind":"strong"');
    expect(rich).toContain('!seq-focus=/refs/x');
    expect(rich).not.toContain('https://x.org');
    expect(rich).not.toContain('‮');
  });

  it('by address: each address a row, located by the caller', async () => {
    const s = await settle(
      spec((b, p) => b.views.tooltip('dbd', { on: `${p}/(range :from 3 :to 6)`, text: 'binds **DNA**' })),
      [SeqRuntime],
    );
    const locate = (a: string): Location | undefined =>
      a === `${P}/(range :from 3 :to 6)` ? rangesLocation(P, [2, 6]) : undefined;
    const ctx = { model: s.model, view: s.view, data: s.store, locate };
    const hit = fresh();
    tooltipsAt(ctx, at(5), hit);
    expect(hit.rows).toEqual([
      {
        label: 'binds DNA',
        value: '',
        rich: [
          { kind: 'text', text: 'binds ' },
          { kind: 'strong', children: [{ kind: 'text', text: 'DNA' }] },
        ],
      },
    ]);
    const miss = fresh();
    tooltipsAt(ctx, at(10), miss);
    expect(miss.rows).toEqual([]);
  });
});

describe('MVS’s schemas', () => {
  it('each row an `elements` form of its schema’s columns alone (a chain file with residue columns is whole chains)', () => {
    const row = { label_asym_id: 'A', beg_label_seq_id: '10', end_label_seq_id: '20', auth_seq_id: '5', tooltip: 'x' };
    const cell = (c: string) => row[c as keyof typeof row];
    expect(rowAddress('/data/c', 'chain', cell)).toBe('/data/c/(elements :label_asym_id A)');
    expect(rowAddress('/data/c', 'residue_range', cell)).toBe(
      '/data/c/(elements :label_asym_id A :beg_label_seq_id 10 :end_label_seq_id 20)',
    );
    expect(rowAddress('/data/c', 'auth_residue', cell)).toBe('/data/c/(elements :auth_seq_id 5)');
    expect(rowAddress('/data/c', 'whole_structure', cell)).toBe('/data/c/(elements)');
    expect(rowAddress('/data/c', 'address', (c) => (c === 'address' ? '(at :label 175)' : undefined))).toBe(
      '/data/c/(at :label 175)',
    );
  });

  it('an MVS annotation’s rows give MVS’s texts (its text_format over field_name), on the rows they locate', async () => {
    // MVS's annotation JSON, as rows: `field_name` note, `text_format` "{}: {score:.1f}"
    const s = await settle(
      spec((b, p) =>
        b.views.tooltip('mvs', {
          target: p,
          schema: 'residue',
          rows: {
            label_asym_id: ['A', 'A'],
            label_seq_id: [5, 9],
            note: ['R175H', 'Y220C'],
            score: [0.912, 0.25],
          } as never,
          field: 'note',
          text: '{}: {score:.1f}',
        }),
      ),
      [SeqRuntime],
    );
    // the caller's structure query, stubbed: residue n at position n − 1
    const locate = (a: string) => {
      const m = /:label_seq_id (\d+)/.exec(a);
      return m ? rangesLocation(P, [Number(m[1]) - 1, Number(m[1])]) : undefined;
    };
    const out = fresh();
    tooltipsAt({ model: s.model, view: s.view, data: s.store, locate }, at(4), out);
    expect(out.rows.map((r) => r.label)).toEqual(['R175H: 0.9']);
    const none = fresh();
    tooltipsAt({ model: s.model, view: s.view, data: s.store, locate }, at(5), none);
    expect(none.rows).toEqual([]);
  });
});
