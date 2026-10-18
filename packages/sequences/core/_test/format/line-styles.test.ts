/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The line style channel in the format (seqspec:format#themes, seqspec:format#layer-types): `ranges`' `line_style` holds
// a line style theme (`uniform` with a `line_style`, or `field_line_style` reading a column of the rows), as `color`
// holds a color theme; it resolves to the layer's themes after color and glyph, before the size channels.
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import {
  compile,
  LineStyles,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
  themeRef,
  validate,
} from '../../src/format/index.js';
import { quiet } from '../../src/format/schema/field-schema.js';
import { KernelLayerTypes } from '../../src/format/tables/kernel.js';
import { files, fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
type B = ReturnType<typeof createSeqViewSpec>;

/** A sequence `s` and ranges `a` with a `line_style` column, and what `more` adds. */
function build(more: (b: B) => void) {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: 's', value: 'MKVLAAGIVGLLLAQ' });
  b.data.annotation({
    ref: 'a',
    geometry: 'ranges',
    target: s,
    value: { start: [0, 3, 6], stop: [2, 5, 9], line_style: ['dashed', 'wavy', 'zigzag'] },
  } as never);
  more(b);
  return b.getDocument() as SeqViewSpecDoc;
}
const codes = (doc: SeqViewSpecDoc) => validate(doc, reg).map((i) => `${i.code} ${i.path.split('/').pop()}`);
const ranges = (layer: Record<string, unknown>, type: 'ranges' | 'points' = 'ranges') =>
  build((b) => {
    const t = b.view({ ref: 'v', label: 'V' }).track({ label: 'T' });
    if (type === 'points') t.points({ data: 'a', ...layer } as never);
    else t.ranges({ ref: 'r', data: 'a', ...layer } as never);
  });
const issuesOf = (ref: unknown) => {
  const ctx = quiet(reg);
  Object.assign(ctx, { lenient: false });
  themeRef('line_style').check(ref, '/t', ctx);
  return ctx.issues.map((i) => `${i.code} ${i.path}`);
};

describe('the line_style channel', () => {
  it('is ranges’ and frame’s, after color and glyph and before the size channels, with no default', () => {
    const styled = KernelLayerTypes.filter((t) => t.channels.some((c) => c.kind === 'line_style')).map((t) => t.name);
    expect(styled).toEqual(['ranges', 'frame']);
    const ch = KernelLayerTypes.find((t) => t.name === 'ranges')?.channels ?? [];
    expect(ch.map((c) => c.name)).toEqual(['color', 'glyph', 'line_style', 'thickness', 'radius']);
    expect(ch[2]).toMatchObject({ kind: 'line_style', unit: null, default: null });
    expect(LineStyles).toEqual(['solid', 'dashed', 'dotted', 'wavy', 'double', 'box', 'highlight']);
  });

  it('holds a constant line style or field_line_style; null is solid', () => {
    expect(codes(ranges({ line_style: { theme: 'uniform', params: { line_style: 'dashed' } } }))).toEqual([]);
    expect(codes(ranges({ line_style: { theme: 'field_line_style' } }))).toEqual([]);
    expect(codes(ranges({ line_style: { theme: 'field_line_style', params: { field: 'kind' } } }))).toEqual([]);
    expect(codes(ranges({ line_style: null }))).toEqual([]);
    for (const s of LineStyles) expect(issuesOf({ theme: 'uniform', params: { line_style: s } })).toEqual([]);
  });

  it('refuses an unknown style, another kind’s theme, an unknown theme, and a type without the channel', () => {
    expect(codes(ranges({ line_style: { theme: 'uniform', params: { line_style: 'zigzag' } } }))).toEqual([
      'svs.line_style.invalid line_style',
    ]);
    expect(issuesOf({ theme: 'uniform', params: { color: 'red' } })).toContain(
      'svs.params.missing /t/params/line_style',
    );
    expect(issuesOf({ theme: 'field_color' })).toEqual(['svs.theme.channel /t/theme']);
    expect(issuesOf({ theme: 'squiggle' })).toEqual(['svs.line_style.invalid /t/theme']);
    expect(codes(ranges({ color: { theme: 'field_line_style' } }))).toEqual(['svs.theme.channel theme']);
    expect(codes(ranges({ line_style: { theme: 'field_line_style' } }, 'points'))).toEqual([
      'svs.params.unknown line_style',
    ]);
  });

  it("resolves to the layer's themes after color and glyph, normalized (field_line_style's field: line_style)", async () => {
    const doc = ranges({
      line_style: { theme: 'field_line_style' },
      thickness: { theme: 'uniform', params: { size: 0.5 } },
    });
    const { set, issues } = await compile(doc, reg, fixtureLoaders(files()));
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    const r = resolveView(doc, set, 'v', reg, supportedAll(reg));
    const L = r.view.layers;
    const k = Array.from(L.ref).indexOf(set.strings.find('r'));
    expect(L.themes[k].map((t) => t?.theme ?? null)).toEqual([
      'field_palette',
      'uniform',
      'field_line_style',
      'uniform',
      null,
    ]);
    expect(L.themes[k][2]).toMatchObject({ params: { field: 'line_style' } });
    // Unknown values in the column are the viewer's to report (svs.render.line_style): the format reads no rows here.
    expect(issues.map((i) => i.code)).not.toContain('svs.render.line_style');
  });
});
