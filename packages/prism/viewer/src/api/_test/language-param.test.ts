/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { compileQuery } from '../../targets/compile.js';
import { parseTarget } from '../../targets/target.js';
import { VmdLanguage } from '../../targets/languages/vmd.js';
import { LanguageParam } from '../target-param.js';

// The `language` param follows the languages the viewer reads (targets.md §1): in this file only VMD (the registry is
// per module, and a test file has its own).

describe('the language param', () => {
  it('is hidden while no language is read, then offers those read, the first by default', async () => {
    const param = LanguageParam();
    expect(param.isHidden).toBe(true);
    expect(param.options).toEqual([]);

    const v = await headlessViewer({ spec: { languages: [VmdLanguage] }, files: [cbs()] });
    expect(param.isHidden).toBe(false);
    expect(param.options).toEqual([['vmd', 'VMD']]);
    expect(PD.getDefaultValues({ language: param })).toEqual({ language: 'vmd' });
    expect(v.api.describe('select.set')?.params.language).toMatchObject({ defaultValue: 'vmd' });

    // Text that names no language is read as the viewer's first.
    expect(compileQuery('resname REA').ok).toBe(true);
    expect(parseTarget('resname REA')).toMatchObject({ ok: true, parsed: { kind: 'query', language: 'vmd' } });
    const r = await v.api.select.set({ target: 'resname REA' });
    expectOk(r);
    expect(await v.api.select.set({ target: 'resn REA', language: 'pymol' as never })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });
});
