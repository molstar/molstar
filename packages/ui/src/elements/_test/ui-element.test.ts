/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { UiElement } from '../ui-element.js';

class TestUiElement extends UiElement<object> {
  protected render(): null {
    return null;
  }
}
customElements.define('msui-test-ui-element', TestUiElement);

afterEach(() => {
  document.body.innerHTML = '';
});

describe('UiElement', () => {
  it('is a root only while it is not inside another UI root, decided at each connect', () => {
    const outer = document.createElement('msui-test-ui-element');
    const inner = document.createElement('msui-test-ui-element');
    document.body.append(outer, inner);
    expect(outer.classList.contains('msui-root')).toBe(true);
    expect(inner.classList.contains('msui-root')).toBe(true);
    // Moved into the other root: no longer a root of its own.
    outer.append(inner);
    expect(inner.classList.contains('msui-root')).toBe(false);
    // Moved out again: a root again.
    document.body.append(inner);
    expect(inner.classList.contains('msui-root')).toBe(true);
  });
});
