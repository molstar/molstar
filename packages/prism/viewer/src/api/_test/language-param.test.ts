/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, languagesOf } from '../../_test/harness.js';
import { compileQuery } from '../../targets/compile.js';
import { parseTarget } from '../../targets/target.js';
import { VmdLanguage } from '../../targets/languages/vmd.js';
import { QueryLanguages } from '../../targets/languages.js';
import { LanguageParam } from '../target-param.js';

// The `language` param follows the languages the viewer reads (targets.md §1): in this file only VMD (the registry is
// per viewer; the page-wide set, the seed of each, is per module and a test file has its own).

describe('the language param', () => {
  it('is hidden while no language is read, then offers those read, the first by default', async () => {
    const set = new QueryLanguages();
    const param = LanguageParam(set);
    expect(param.isHidden).toBe(true);
    expect(param.options).toEqual([]);

    const v = await headlessViewer({ spec: { languages: [VmdLanguage] }, files: [cbs()] });
    set.register([VmdLanguage]);
    expect(param.isHidden).toBe(false);
    expect(param.options).toEqual([['vmd', 'VMD']]);
    expect(PD.getDefaultValues({ language: param })).toEqual({ language: 'vmd' });
    expect(v.api.describe('select.set')?.params.language).toMatchObject({ defaultValue: 'vmd' });

    // Text that names no language is read as the viewer's first.
    const languages = languagesOf(v);
    expect(compileQuery('resname REA', undefined, languages).ok).toBe(true);
    expect(parseTarget('resname REA', { languages })).toMatchObject({
      ok: true,
      parsed: { kind: 'query', language: 'vmd' },
    });
    const r = await v.api.select.set({ target: 'resname REA' });
    expectOk(r);
    expect(await v.api.select.set({ target: 'resn REA', language: 'pymol' as never })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });
});
