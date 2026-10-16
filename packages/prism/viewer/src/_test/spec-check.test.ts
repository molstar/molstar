/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { checkSpec, isPartId } from '../spec-check.js';

// The spec check (viewer-spec.md §6): part ids are told from command ids by their kind.

describe('the spec check', () => {
  it('knows every kind of part, a full-app page too', () => {
    for (const kind of ['tab', 'section', 'row', 'toolbar', 'overlay', 'page', 'tool', 'skill', 'settings']) {
      expect(isPartId(`ext.${kind}.name`), kind).toBe(true);
    }
    expect(isPartId('ext.nothing.name')).toBe(false);
  });

  it('checks the owner of an overridden page like that of any other part', () => {
    const spec = { ...DefaultPrismViewerSpec(), overrides: { 'nobody.page.inspector': null } };
    expect(() => checkSpec(spec)).toThrow("override 'nobody.page.inspector': there is no extension 'nobody'");
  });
});
