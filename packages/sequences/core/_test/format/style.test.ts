/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's `style` as a preset name or an inline preset (seqspec:format#sections): the validator's warnings, each
// naming what is allowed, and resolve's mapping of the format's keys onto the viewer's render props.

import { describe, expect, it } from 'vitest';
import {
  compile,
  resolveStyle,
  resolveView,
  SectionStyleKeys,
  supportedAll,
  validate,
} from '../../src/format/index.js';
import { type J, kernel, minimalDoc } from './fixtures.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

const withStyle = (style: unknown): J => {
  const d = minimalDoc();
  d.root.children[1].children[0].params = { style };
  return d;
};
const PATH = '/root/children/1/children/0/params/style';
const issues = (style: unknown, version = '1.0.0') => {
  const d = withStyle(style);
  d.version = version;
  return validate(d, kernel()).map((i) => [i.code, i.severity, i.path.slice(PATH.length), i.message]);
};

describe("validating a section's style", () => {
  it('a preset name, as in 1.0, and an inline preset with each key valid: no issues', () => {
    expect(issues('msa')).toEqual([]);
    expect(issues('org.example:dense')).toEqual([]);
    const each: Record<string, unknown> = {
      preset: 'msa',
      column_heads: 'own_row',
      data_floor: 0,
      row_height: 18,
      header_width: 0,
      vertical_nav: 'both',
      wrapped_group_labels: 'inline_start',
      row_borders: true,
      hover_emphasis: false,
    };
    expect(Object.keys(each).slice(1)).toEqual(Object.keys(SectionStyleKeys));
    for (const [k, v] of Object.entries(each)) expect(issues({ [k]: v }), k).toEqual([]);
    expect(issues(each)).toEqual([]);
    expect(issues({})).toEqual([]);
    expect(
      issues({ preset: null, column_heads: 'shared', vertical_nav: 'none', wrapped_group_labels: 'none' }),
    ).toEqual([]);
  });

  it('an unknown key is a warning naming every key; the rest is kept', () => {
    expect(issues({ preset: 'msa', colum_heads: 'own_row', data_floor: 0 })).toEqual([
      [
        'svs.section.style',
        'warning',
        '/colum_heads',
        'unknown style key "colum_heads"; the keys: preset, column_heads, data_floor, row_height, header_width, ' +
          'vertical_nav, wrapped_group_labels, row_borders, hover_emphasis; it is ignored, the rest of the style applies',
      ],
    ]);
    // A newer minor's key is the version's warning, as any unknown param there (seqspec:format#versioning).
    expect(issues({ ring_gap: 4 }, '1.1.0').map((i) => i.slice(0, 3))).toEqual([
      ['svs.version.newer_minor', 'warning', ''],
      ['svs.version.newer_minor', 'warning', '/ring_gap'],
    ]);
  });

  it('a bad value is a warning naming the values allowed', () => {
    const bad = (key: string, value: unknown) => issues({ [key]: value }).map((i) => [i[2], i[3]]);
    const tail = '; it is ignored, the rest of the style applies';
    expect(bad('column_heads', 'ownRow')).toEqual([
      ['/column_heads', `column_heads "ownRow" is not one of "shared", "own_row"${tail}`],
    ]);
    expect(bad('data_floor', -1)).toEqual([['/data_floor', `data_floor -1 is not a number of px, at least 0${tail}`]]);
    expect(bad('row_height', 0)).toEqual([['/row_height', `row_height 0 is not a number of px, at least 1${tail}`]]);
    expect(bad('header_width', '150px')).toEqual([
      ['/header_width', `header_width "150px" is not a number of px, at least 0${tail}`],
    ]);
    expect(bad('vertical_nav', 'map')).toEqual([
      ['/vertical_nav', `vertical_nav "map" is not one of "scrollbar", "overview", "both", "none"${tail}`],
    ]);
    expect(bad('wrapped_group_labels', 'aboveEnd')).toEqual([
      [
        '/wrapped_group_labels',
        `wrapped_group_labels "aboveEnd" is not one of "above_end", "inline_start", "none"${tail}`,
      ],
    ]);
    expect(bad('row_borders', 'yes')).toEqual([['/row_borders', `row_borders "yes" is not true or false${tail}`]]);
    expect(bad('preset', 3)).toEqual([['/preset', `preset 3 is not a preset name or null${tail}`]]);
  });

  it('a style neither a name nor an object is an error, as a non-string was in 1.0', () => {
    expect(issues(3).map((i) => i.slice(0, 2))).toEqual([['svs.params.type', 'error']]);
    expect(issues(['msa']).map((i) => i.slice(0, 2))).toEqual([['svs.params.type', 'error']]);
  });
});

describe("resolving a section's style", () => {
  it("maps the format's keys onto the render props, the valid ones only", () => {
    expect(resolveStyle('msa')).toEqual({ preset: 'msa', render: null });
    expect(resolveStyle(null)).toEqual({ preset: null, render: null });
    expect(resolveStyle({ data_floor: 0 })).toEqual({ preset: null, render: { dataFloor: 0 } });
    expect(
      resolveStyle({
        preset: 'msa',
        column_heads: 'own_row',
        data_floor: 32,
        row_height: 16,
        header_width: 120,
        vertical_nav: 'overview',
        wrapped_group_labels: 'above_end',
        row_borders: false,
        colum_heads: 'own_row',
        row_borders_: true,
      }),
    ).toEqual({
      preset: 'msa',
      render: {
        dataColumnHeads: 'ownRow',
        dataFloor: 32,
        rowHeight: 16,
        headerWidth: 120,
        verticalNav: 'overview',
        wrappedGroupLabels: 'aboveEnd',
        rowBorders: false,
      },
    });
    expect(resolveStyle({ preset: 3, column_heads: 'ownRow', data_floor: -1 })).toEqual({
      preset: null,
      render: null,
    });
  });

  it("the resolved view carries the preset's name and the inline keys per section", async () => {
    const reg = registries();
    for (const [style, preset, render] of [
      ['strip', 'strip', null],
      [{ preset: 'msa', column_heads: 'own_row', data_floor: 0 }, 'msa', { dataColumnHeads: 'ownRow', dataFloor: 0 }],
      [{ inline_start: true, row_height: 30 }, null, { rowHeight: 30 }],
    ] as const) {
      const doc = withStyle(style);
      const { set, issues } = await compile(doc, reg, fixtureLoaders());
      expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
      const { view } = resolveView(doc, set, 'v', reg, supportedAll(reg));
      expect(view.sections.style[0] < 0 ? null : set.strings.get(view.sections.style[0])).toBe(preset);
      expect(view.sections.render[0]).toEqual(render);
    }
  });
});
