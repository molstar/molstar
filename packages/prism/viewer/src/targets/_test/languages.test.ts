/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, languagesOf } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { PrismViewer } from '../../viewer.js';
import { compileQuery } from '../compile.js';
import { queryLanguages, registerQueryLanguages } from '../languages.js';
import { PyMolLanguage } from '../languages/pymol.js';
import { JmolLanguage } from '../languages/jmol.js';
import { VmdLanguage } from '../languages/vmd.js';
import { parseTarget } from '../target.js';

// The query languages a viewer reads are its spec's (targets/languages.ts): with none, keywords and address lists still
// read and text says so; with PyMOL only, another language is refused by name. In order: the page-wide set is per module; each viewer has its own on top
// of it (README decision 29).

describe('query languages', () => {
  it('reads keywords and address lists with none, and says text needs a language', () => {
    expect(queryLanguages()).toEqual([]);
    expect(parseTarget('ligand')).toMatchObject({ ok: true, parsed: { kind: 'keyword' } });
    expect(parseTarget('A:120-130')).toMatchObject({ ok: true, parsed: { kind: 'addresses' } });
    expect(compileQuery('chain A and resn HIS')).toMatchObject({
      ok: false,
      message: 'It is not a keyword or an address list, and this viewer reads no query language.',
    });
  });

  it('reads the languages registered, and refuses another by name', () => {
    registerQueryLanguages([PyMolLanguage]);
    registerQueryLanguages([PyMolLanguage]);
    expect(queryLanguages().map((l) => l.id)).toEqual(['pymol']);
    expect(compileQuery('chain A and resn HIS').ok).toBe(true);
    const vmd = compileQuery('resname HIS', 'vmd');
    expect(vmd).toMatchObject({ ok: false, message: 'This viewer does not read VMD.' });
    expect(!vmd.ok && vmd.hint).toMatch(/^Write it in PyMOL, or as a keyword or an address list/);
  });

  it("gives each viewer its own set: the page-wide defaults, then its spec's", async () => {
    const a = await createPrismViewerHeadless({ spec: { ...DefaultPrismViewerSpec(), languages: [VmdLanguage] } });
    const b = await createPrismViewerHeadless({ spec: { ...DefaultPrismViewerSpec(), languages: [JmolLanguage] } });
    try {
      const ids = (v: PrismViewer) => languagesOf(v).list.map((l) => l.id);
      expect(ids(a)).toEqual(['pymol', 'vmd']);
      expect(ids(b)).toEqual(['pymol', 'jmol']);
      // Nothing leaks to the page-wide set, nor to the other viewer.
      expect(queryLanguages().map((l) => l.id)).toEqual(['pymol']);
      expect(languagesOf(a).get('jmol')).toBeUndefined();
      expect(languagesOf(b).get('vmd')).toBeUndefined();

      // What they read and offer follows: VMD text on A only, Jmol text on B only.
      expect(compileQuery('resname HIS', 'vmd', languagesOf(a)).ok).toBe(true);
      expect(compileQuery('resname HIS', 'vmd', languagesOf(b))).toMatchObject({
        ok: false,
        message: 'This viewer does not read VMD.',
      });
      await a.api.load.files({ files: [cbs()] });
      await b.api.load.files({ files: [cbs()] });
      expect(await a.api.select.set({ target: 'resname REA', language: 'vmd' })).toMatchObject({ ok: true });
      expect(await b.api.select.set({ target: 'resname REA', language: 'vmd' as never })).toMatchObject({
        ok: false,
        code: 'invalid',
      });
      expect(a.api.describe('select.set')?.params.language).toMatchObject({
        options: [
          ['pymol', 'PyMOL'],
          ['vmd', 'VMD'],
        ],
      });
      expect(b.api.describe('select.set')?.params.language).toMatchObject({
        options: [
          ['pymol', 'PyMOL'],
          ['jmol', 'Jmol'],
        ],
      });
    } finally {
      await a.dispose();
      await b.dispose();
    }
  });
});
