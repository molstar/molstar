/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom
import { beforeEach, describe, expect, it } from 'vitest';
import { assistantPrefsPrefix, localPrefs } from '../settings-model.js';

// The assistant's choices (model, own endpoints, allowed domains) are kept per viewer, by its persist key; the app's
// key keeps the original prefix, no key is `prism-default` (viewer README decision 29).

describe('the assistant preferences of several viewers on a page', () => {
  beforeEach(() => localStorage.clear());

  it("keep the original prefix for the app's key, prism-default's for no key, and a prefix of their own for another", () => {
    expect(assistantPrefsPrefix('prism-viewer')).toBe('prism-assistant.');
    // No key is `prism-default` (as the kept records are): not the app's choices.
    expect(assistantPrefsPrefix()).toBe('prism-assistant@prism-default.');
    expect(assistantPrefsPrefix('')).toBe('prism-assistant@prism-default.');
    expect(assistantPrefsPrefix('left')).toBe('prism-assistant@left.');
    expect(assistantPrefsPrefix('left')).not.toBe(assistantPrefsPrefix('right'));
  });

  it('do not see each other', () => {
    const left = localPrefs(assistantPrefsPrefix('left'));
    const right = localPrefs(assistantPrefsPrefix('right'));
    left.set('domains', '["a.org"]');
    expect(left.get('domains')).toBe('["a.org"]');
    expect(right.get('domains')).toBeUndefined();
    // Never under the keys' prefix (the person's keys are shared, `prism-assistant.key.<provider>`).
    expect([...Object.keys(localStorage)].some((k) => k.startsWith('prism-assistant.key.'))).toBe(false);
  });
});
