/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The conformance suite over every built-in scheme and both token tables (_spec/README.md "Themes", Schemes; "Render",
// Tokens): ink at 4.5:1 on every cell color; missing apart from every ramp's low end; the category tokens and the safe
// schemes distinct under the color-vision check; a ramp monotone in lightness, adjacent stops 0.03 apart; a color in
// each table for each scheme (its dark variant); the centromere at 3:1 against the surface.

import { checkColorVision } from '@molstar/core/util/color/cvd';
import { Oklab } from '@molstar/core/util/color/spaces/oklab';
import { DarkTokens, inkOn, LightTokens, type SeqTokens } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { schemeProvenance } from '../provenance.js';
import { BuiltInSchemes, compileChannel, contrastRgba, namedList, type SchemeEntry, toRgb } from '../index.js';
import { builtInSchemes } from '../schemes.js';

const Tables: readonly SeqTokens[] = [LightTokens, DarkTokens];
const registry = builtInSchemes();
// what each subject is drawn over: letters, numbers, categories, the stain of a band
const datums = (e: SchemeEntry) => {
  const out: { letter: string; value: number; field: (n: string) => string | number | undefined }[] = [];
  const add = (letter: string, value: number, v?: string) =>
    out.push({ letter, value, field: (n) => (n === 'value' ? value : v) });
  if (e.subject === 'protein') for (const l of 'ACDEFGHIKLMNPQRSTVWY-') add(l, Number.NaN);
  else if (e.subject === 'nucleotide') for (const l of 'ACGTUN-') add(l, Number.NaN);
  else if (e.subject === 'numbers') for (let i = 0; i <= 10; i++) add('', i * 10);
  else if (e.subject === 'secondary') for (const k of ['helix', 'strand', 'turn']) add('', Number.NaN, k);
  else if (e.subject === 'bands')
    for (const k of ['gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100', 'acen', 'gvar', 'stalk']) add('', Number.NaN, k);
  else for (const k of ['a', 'b', 'c', 'd', 'e']) add('', Number.NaN, k);
  return out;
};

/** The colors a scheme draws in a table (measures unregistered: column schemes are drawn as their default there). */
function drawn(e: SchemeEntry, tokens: SeqTokens): number[] {
  const alphabet = e.subject === 'protein' ? 'protein' : e.subject === 'nucleotide' ? 'dna' : undefined;
  const measures = new Map((e.needs ?? []).map((m) => [m, { name: m, read: () => 0.5 }]));
  const ch = compileChannel('color', { kind: 'scheme', params: { name: e.name } }, 'token:missing', {
    tokens,
    schemes: registry,
    measures,
    ...(alphabet ? { alphabet } : {}),
  });
  const out: number[] = [];
  for (const d of datums(e)) {
    const c = ch.at({ ...d, start: 0, end: 1 });
    if (typeof c === 'number') out.push(c);
  }
  return out;
}

describe('the built-in schemes in both tables', () => {
  for (const e of BuiltInSchemes)
    it(`${e.name}: colors in each table, ink at 4.5:1 on each, a provenance`, () => {
      expect(schemeProvenance(e.name)?.length ?? 0).toBeGreaterThan(10);
      for (const tokens of Tables) {
        const colors = drawn(e, tokens);
        expect(colors.length, `${e.name} draws in ${tokens.theme}`).toBeGreaterThan(0);
        for (const c of colors)
          expect(contrastRgba(inkOn(c, tokens), c), `${e.name} ${tokens.theme}`).toBeGreaterThanOrEqual(4.5);
      }
    });

  it('the safe categorical and letter schemes are distinct under the color-vision check, in both tables', () => {
    for (const e of BuiltInSchemes.filter(
      (x) => x.cvdSafe && (x.subject === 'categories' || x.subject === 'nucleotide'),
    ))
      for (const tokens of Tables) {
        // the scheme's own colors: a letter it does not list takes the missing token
        const colors = [...new Set(drawn(e, tokens))].filter((c) => c !== tokens.colors['seq-missing']).map(toRgb);
        expect(checkColorVision(colors, 0.06).ok, `${e.name} ${tokens.theme}`).toBe(true);
      }
  });

  it('a scheme that is not safe declares a safe substitute that fits its subject', () => {
    for (const e of BuiltInSchemes.filter((x) => !x.cvdSafe && x.substitute)) {
      const s = registry.get(e.substitute!)!;
      expect(s.subject).toBe(e.subject);
    }
  });
});

describe('the token tables', () => {
  for (const tokens of Tables)
    it(`${tokens.theme}: categories pass the check; ramps are monotone; missing stands apart; the centromere reads`, () => {
      const cats = namedList('categories', tokens)!.map(toRgb);
      const check = checkColorVision(cats, 0.06);
      expect(check.ok, JSON.stringify(check.worst)).toBe(true);
      const ramp = namedList('ramp', tokens)!;
      const ls = ramp.map((c) => Oklab.lightness(toRgb(c)));
      const dir = Math.sign(ls[ls.length - 1] - ls[0]);
      for (let i = 1; i < ls.length; i++) expect(dir * (ls[i] - ls[i - 1])).toBeGreaterThanOrEqual(0.03);
      const missing = toRgb(tokens.colors['seq-missing']);
      for (const name of ['ramp', 'diverging', 'viridis', 'cividis', 'magma', 'blues', 'rd-bu']) {
        const low = toRgb(namedList(name, tokens)![0]);
        expect(Oklab.colorDistance(low, missing), `${name}’s low end`).toBeGreaterThanOrEqual(0.05);
      }
      expect(contrastRgba(tokens.colors['seq-centromere'], tokens.colors['msui-bg'])).toBeGreaterThanOrEqual(3);
      expect(contrastRgba(tokens.colors['seq-centromere'], tokens.colors['seq-track'])).toBeGreaterThanOrEqual(3);
    });
});
