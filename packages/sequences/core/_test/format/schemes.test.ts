/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Color schemes (seqspec:format#color-schemes): each v1 scheme against its published definition for a few letters or
// stops, the dark variants, the inline scheme's validation, what fits which subject, and resolve drawing a scheme of
// another subject as the layer's default.

import { KernelSchemeFits } from '@molstar/seq-core/builder';
import type { Issue, SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  applyPatch,
  type ChannelKind,
  Color,
  ColorLists,
  compile,
  KernelColorSchemes,
  LetterDictionaries,
  type Patch,
  resolveView,
  type SchemeParams,
  type SchemeSubject,
  schemeLegend,
  schemesFor,
  schemeTheme,
  supportedAll,
  type ThemeRef,
  themeRef,
} from '../../src/format/index.js';
import { darkRamp } from '../../src/format/palettes/schemes.js';
import { quiet } from '../../src/format/schema/field-schema.js';
import { CASES, case3 } from './cases.js';
import { files, fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
const theme = (p: SchemeParams, dark = false) => schemeTheme(reg, p, { dark }) as ThemeRef;
const letters = (name: string, dark = false) =>
  (theme({ name }, dark).params as { colors: Record<string, string> }).colors;
const palette = (name: string, dark = false) =>
  (theme({ name }, dark).params as { palette: { colors: readonly string[] } }).palette.colors;
const issuesOf = (ref: unknown, kind: ChannelKind = 'color') => {
  const ctx = quiet(reg);
  Object.assign(ctx, { lenient: false });
  themeRef(kind).check(ref, '/t', ctx);
  return ctx.issues;
};
const scheme = (params: Record<string, unknown>) => ({ theme: 'scheme', params });
const patched = (doc: ReturnType<typeof case3>, ops: Patch['ops']) =>
  applyPatch(doc, { format: 'seq_view_spec_patch', ops }, reg).doc;

describe('the v1 schemes against their published definitions', () => {
  it('17, each named once, with a label, a description and its provenance', () => {
    const names = KernelColorSchemes.map((e) => e.name);
    expect(names).toHaveLength(17);
    expect(new Set(names).size).toBe(17);
    for (const e of KernelColorSchemes) for (const s of [e.label, e.description, e.provenance]) expect(s).not.toBe('');
  });

  it("protein letters: Zappo, Taylor, hydrophobicity (Jalview's colors)", () => {
    const z = letters('zappo');
    expect([z.A, z.I, z.K, z.D, z.S, z.G, z.C, z.W]).toEqual([
      '#ffafaf',
      '#ffafaf',
      '#6464ff',
      '#ff0000',
      '#00ff00',
      '#ff00ff',
      '#ffff00',
      '#ffc800',
    ]);
    // Zappo's groups give exactly Jalview's dictionary.
    expect(z).toEqual(
      Object.fromEntries(Object.entries(LetterDictionaries.Zappo).map(([k, c]) => [k, Color.toHex(c)])),
    );
    const t = letters('taylor');
    expect([t.A, t.R, t.D, t.W]).toEqual(['#ccff00', '#0000ff', '#ff0000', '#00ccff']);
    const h = letters('hydrophobicity');
    expect([h.I, h.R, h.K]).toEqual(['#ff0000', '#0000ff', '#0000ff']); // Kyte–Doolittle's ends: I most hydrophobic
  });

  it("nucleotides: Jalview's nucleotide and purine/pyrimidine, and Okabe–Ito's", () => {
    const n = letters('nucleotide');
    expect([n.A, n.C, n.G, n.T, n.U]).toEqual(['#64f73f', '#ffb340', '#eb413c', '#3c88ee', '#3c88ee']);
    const p = letters('purine_pyrimidine');
    expect([p.A, p.G, p.C, p.T]).toEqual(['#ff83fa', '#ff83fa', '#40e0d0', '#40e0d0']);
    const cb = letters('nucleotide_cb');
    expect([cb.A, cb.C, cb.G, cb.T]).toEqual(['#009e73', '#0072b2', '#f0e442', '#d55e00']);
  });

  it("Clustal X is the kernel's table; percent identity Jalview's thresholds over the consensus", () => {
    expect(theme({ name: 'clustal_x' })).toEqual({ theme: 'clustal_x' });
    expect(theme({ name: 'percent_identity' })).toEqual({ theme: 'percent_identity' });
    const table = reg.get('themes', 'percent_identity')?.expand?.({}, 'color') as ThemeRef;
    const rules = (table.params as { rules: { when: { measure: unknown; min?: number }[]; color: string }[] }).rules;
    expect(rules.map((r) => r.color)).toEqual(['#6464ff', '#9999ff', '#ccccff']);
    expect(rules.map((r) => Math.round((r.when[1].min ?? 0) * 100))).toEqual([80, 60, 40]);
    expect(rules[0].when[0].measure).toEqual({ rank: 'letter' });
  });

  it("values and categories: MVS's lists, Okabe–Ito's seven colors", () => {
    expect(palette('viridis')).toEqual(ColorLists.Viridis.map(Color.toHex));
    expect(palette('rd_bu')[5]).toBe('#f7f7f7');
    expect(palette('okabe_ito')).toEqual(['#e69f00', '#56b4e9', '#009e73', '#f0e442', '#0072b2', '#d55e00', '#cc79a7']);
    expect(theme({ name: 'viridis', field: 'plddt', domain: [0, 100], reverse: true }).params).toMatchObject({
      field: 'plddt',
      palette: { kind: 'continuous', value_domain: [0, 100], reverse: true },
    });
  });

  it("legends: a letters scheme's labeled groups, Clustal X's groups, overrides recolor them", () => {
    expect(schemeLegend(reg, { name: 'purine_pyrimidine' }, false)).toEqual([
      { label: 'purine', color: '#ff83fa' },
      { label: 'pyrimidine', color: '#40e0d0' },
    ]);
    expect(schemeLegend(reg, { name: 'clustal_x' }, false).map((x) => x.label)).toContain('hydrophobic');
    expect(schemeLegend(reg, { name: 'zappo', overrides: { C: '#f0c' } }, false).at(-1)).toEqual({
      label: 'cysteine',
      color: '#f0c',
    });
  });

  it("legends: a theme scheme's dark variant never lists the light theme's colors", () => {
    // ESPript's and percent identity's rules are labeled: their theme lists them, in the colors drawn on either skin.
    for (const name of ['espript', 'percent_identity']) {
      expect(schemeLegend(reg, { name }, false)).toEqual([]);
      expect(schemeLegend(reg, { name }, true)).toEqual([]);
    }
    // Every theme scheme with a dark variant: its dark legend is `darkLegend`, never `legend`.
    for (const e of KernelColorSchemes) {
      const def = e.scheme;
      if (def.kind !== 'theme' || !def.dark) continue;
      expect(schemeLegend(reg, { name: e.name }, true)).toEqual([...(def.darkLegend ?? [])]);
    }
    const custom = {
      kind: 'theme',
      alphabets: ['protein'],
      columns: true,
      theme: { theme: 'espript' },
      dark: { theme: 'espript', params: { identical: '#f0f0f0' } },
      legend: [{ label: 'x', color: '#000000' }],
      darkLegend: [{ label: 'x', color: '#f0f0f0' }],
    } as const;
    expect(schemeLegend(reg, { custom }, false)).toEqual([{ label: 'x', color: '#000000' }]);
    expect(schemeLegend(reg, { custom }, true)).toEqual([{ label: 'x', color: '#f0f0f0' }]);
  });

  it('legends: a letters scheme with no labeled groups lists its letters, those of one color together', () => {
    const taylor = schemeLegend(reg, { name: 'taylor' }, false);
    expect(taylor.length).toBeGreaterThanOrEqual(20);
    expect(taylor.some((x) => x.label.split(' ').includes('W'))).toBe(true);
    const hydro = schemeLegend(reg, { name: 'hydrophobicity' }, false);
    expect(hydro.length).toBeGreaterThan(0);
    // Every letter once, across the items.
    const all = hydro.flatMap((x) => x.label.split(' '));
    expect(new Set(all).size).toBe(all.length);
    expect(all).toContain('I');
  });
});

describe('dark variants', () => {
  it('a light-low sequential ramp is reversed; a dark-low one stays', () => {
    expect(palette('blues', true)).toEqual([...palette('blues')].reverse());
    expect(palette('viridis', true)).toEqual(palette('viridis'));
  });

  it('a diverging ramp: each half reversed, the center kept, each hue on its side', () => {
    const light = palette('rd_bu');
    const dark = palette('rd_bu', true);
    expect(dark[0]).toBe(light[4]);
    expect(dark[4]).toBe(light[0]);
    expect(dark[5]).toBe(light[5]);
    expect(dark[10]).toBe(light[6]);
    expect(darkRamp('diverging', ['#f00', '#fcc', '#ccf', '#00f'])).toEqual(['#fcc', '#f00', '#00f', '#ccf']);
  });

  it('letters and categories keep their colors; percent identity has its own; a custom dark variant applies', () => {
    expect(letters('zappo', true)).toEqual(letters('zappo'));
    expect(palette('tableau10', true)).toEqual(palette('tableau10'));
    const dark = theme({ name: 'percent_identity' }, true);
    expect((dark.params as { rules: { color: string }[] }).rules[0].color).toBe('#8080ff');
    const custom = {
      kind: 'letters' as const,
      alphabets: ['dna' as const],
      colors: { A: '#0a0' },
      dark: { colors: { A: '#5f5' } },
    };
    expect((theme({ custom }, true).params as { colors: Record<string, string> }).colors.A).toBe('#5f5');
    expect(
      (theme({ custom: { ...custom, dark: 'same' } }, true).params as { colors: Record<string, string> }).colors.A,
    ).toBe('#0a0');
  });
});

describe('a custom scheme', () => {
  it("groups, a default for the alphabet's other letters, overrides last", () => {
    const custom = {
      kind: 'letters' as const,
      alphabets: ['dna' as const],
      groups: [{ letters: 'AG', color: '#f00', label: 'purine' }],
      default: '#ccc',
    };
    const c = (theme({ custom, overrides: { G: '#00f' } }).params as { colors: Record<string, string> }).colors;
    expect([c.A, c.G, c.C, c.T, c.N, c.R]).toEqual(['#f00', '#00f', '#ccc', '#ccc', '#ccc', '#ccc']);
    expect(c.U).toBeUndefined(); // not a DNA letter
    expect(issuesOf(scheme({ custom, overrides: { G: '#00f' } }))).toEqual([]);
  });

  it('instructive errors: a letter the alphabet lacks, a bad color, a ramp of one, overrides on values', () => {
    const bad = issuesOf(scheme({ custom: { kind: 'letters', alphabets: ['dna'], colors: { A: '#0a0', J: '#00f' } } }));
    expect(bad.map((i) => `${i.code} ${i.path}`)).toEqual(['svs.scheme.letter /t/params/custom/colors/J']);
    expect(bad[0].message).toContain('dna (ABCDGHKMNRSTVWY)');
    const color = issuesOf(scheme({ custom: { kind: 'sequential', colors: ['#fff', 'bluish'] } }));
    expect(color.map((i) => `${i.code} ${i.path}`)).toEqual(['svs.color.invalid /t/params/custom/colors/1']);
    expect(issuesOf(scheme({ custom: { kind: 'diverging', colors: ['#fff'] } }))[0].message).toContain('at least two');
    expect(issuesOf(scheme({ custom: { kind: 'letters', alphabets: ['amino'] } }))[0].message).toContain(
      'protein, dna',
    );
    expect(issuesOf(scheme({ name: 'viridis', overrides: { A: '#f00' } }))[0].message).toContain('custom scheme');
    expect(issuesOf(scheme({ name: 'zappo', overrides: { '-': '#f00' } })).map((i) => i.code)).toEqual([
      'svs.scheme.letter',
    ]);
    expect(issuesOf(scheme({ name: 'zappo', custom: { kind: 'sequential', colors: ['#fff', '#000'] } }))[0].code).toBe(
      'svs.params.exclusive',
    );
    expect(issuesOf(scheme({}))[0].code).toBe('svs.params.missing');
  });

  it("an unknown name lists the registered ones; an extension's unknown here is a warning", () => {
    const [unknown] = issuesOf(scheme({ name: 'rainbow' }));
    expect(unknown.code).toBe('svs.scheme.unknown');
    expect(unknown.message).toContain('clustal_x, zappo');
    const [ext] = issuesOf(scheme({ name: 'org.example:brand' }));
    expect([ext.code, ext.severity]).toEqual(['svs.render.theme', 'warning']);
  });
});

describe("the builder's twin", () => {
  it("names the kernel's schemes in order, each fitting what it fits here", () => {
    const fits = KernelColorSchemes.map((e) => {
      const d = e.scheme;
      if (d.kind === 'letters') return [e.name, { letters: d.alphabets }];
      if (d.kind === 'theme')
        return [e.name, d.columns ? { letters: d.alphabets, columns: true } : { letters: d.alphabets }];
      return [e.name, d.kind === 'categorical' ? 'categories' : 'values'];
    });
    expect(Object.entries(KernelSchemeFits)).toEqual(fits);
  });
});

describe('what fits', () => {
  const names = (s: SchemeSubject) => schemesFor(reg, s).map((e) => e.name);
  it('letters by alphabet, column tables only on held columns; values; categories', () => {
    const protein = ['clustal_x', 'zappo', 'taylor', 'hydrophobicity'];
    expect(names({ kind: 'letters', alphabet: 'protein', columns: false })).toEqual(protein);
    expect(names({ kind: 'letters', alphabet: 'protein', columns: true })).toEqual([
      ...protein,
      'percent_identity',
      'espript',
    ]);
    expect(names({ kind: 'letters', alphabet: 'dna', columns: false })).toEqual([
      'nucleotide',
      'purine_pyrimidine',
      'nucleotide_cb',
    ]);
    expect(names({ kind: 'values' })).toEqual(['viridis', 'cividis', 'magma', 'blues', 'rd_bu']);
    expect(names({ kind: 'categories' })).toEqual(['category10', 'okabe_ito', 'tableau10']);
  });
});

describe('resolve', () => {
  const all = supportedAll(reg);
  const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
  async function color(doc: ReturnType<typeof case3>, view: string, layer: string) {
    const { set, issues } = await compile(doc, reg, fixtureLoaders(files()));
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    const r = resolveView(doc, set, view, reg, all);
    const l = Array.from(r.view.layers.ref).indexOf((set as SeqSet).strings.find(layer));
    return { theme: r.view.layers.themes[l][0], issues: r.issues };
  }

  it("a scheme that fits is kept for the viewer; one of another subject is the layer's default, reported", async () => {
    const set = (c: unknown): Patch['ops'] => [{ op: 'set', ref: 'family.rows', params: { color: c } }];
    const ok = await color(patched(case3(), set(scheme({ name: 'zappo' }))), 'family', 'family.rows');
    expect(ok.theme?.theme).toBe('scheme');
    expect(codes(ok.issues)).toEqual([]);
    const off = await color(patched(case3(), set(scheme({ name: 'viridis' }))), 'family', 'family.rows');
    expect(off.theme).toMatchObject({ theme: 'letter' });
    expect(codes(off.issues)).toEqual(['svs.render.theme']);
    // Percent identity reads only ranks and shares: a nucleotide alignment's held columns too.
    const rna = patched(CASES.case4(), [
      { op: 'set', ref: 'rna.rows', params: { color: scheme({ name: 'percent_identity' }) } },
    ]);
    const pi = await color(rna, 'rna', 'rna.rows');
    expect(pi.theme?.theme).toBe('scheme');
    expect(codes(pi.issues)).toEqual([]);
  });
});
