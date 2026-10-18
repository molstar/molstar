/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';

describe('script languages', () => {
  it('only mol-script is available until a transpiler module is imported', async () => {
    vi.resetModules();
    const { Script } = await import('@molstar/model/script/script');
    expect(Script.getAvailableLanguages()).toEqual(['mol-script']);
    expect(() => Script.toExpression({ language: 'pymol', expression: 'resn ALA' })).toThrow(
      "Script language 'pymol' is not available in this build",
    );
    expect(() => Script.toExpression({ language: 'vmd', expression: 'resname ALA' })).toThrow(
      "Script language 'vmd' is not available in this build",
    );
    expect(() => Script.toExpression({ language: 'jmol', expression: 'ALA' })).toThrow(
      "Script language 'jmol' is not available in this build",
    );
  });

  it('importing a transpiler module enables its language only', async () => {
    vi.resetModules();
    const { Script } = await import('@molstar/model/script/script');
    await import('@molstar/query-language/transpilers/pymol');
    expect(Script.getAvailableLanguages()).toEqual(['mol-script', 'pymol']);
    expect(() => Script.toExpression({ language: 'pymol', expression: 'resn ALA' })).not.toThrow();
    expect(() => Script.toExpression({ language: 'vmd', expression: 'resname ALA' })).toThrow(/not available/);
  });

  it('transpilers/all enables every language', async () => {
    vi.resetModules();
    const { Script } = await import('@molstar/model/script/script');
    await import('@molstar/query-language/transpilers/all');
    expect(Script.getAvailableLanguages().sort()).toEqual(['jmol', 'mol-script', 'pymol', 'vmd']);
  });
});
