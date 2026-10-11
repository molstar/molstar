/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from '../../components/button/_generated/button.module.css.js';
import { defineElements, MspElement } from '../../dom/element.js';
import { defineStyles, styleText } from '../registry.js';

class StyledElement extends MspElement<string> {
  static tag = 'msp-test-styled';
  protected render(text: string) {
    return <span class={classes.button}>{text}</span>;
  }
}
defineElements([StyledElement]);

const tick = () => new Promise<void>((resolve) => queueMicrotask(resolve));

describe('styles', () => {
  it('compiles modules to prefixed class names in the components layer', () => {
    expect(classes.button).toBe('msp-button-button');
    expect(classes.iconOnly).toBe('msp-button-iconOnly');
    const text = styleText();
    expect(text.startsWith('@layer msp.reset, msp.theme, msp.components, msp.overrides;')).toBe(true);
    expect(text).toContain('@layer msp.components {\n.msp-button-button');
    expect(text).toContain('--msp-color-bg: #eeece7');
  });

  it('adopts the registered styles when an element connects, and follows later registrations', async () => {
    const el = document.createElement('msp-test-styled') as StyledElement;
    el.model = 'hi';
    document.body.append(el);
    const style = document.head.querySelector('style[data-msp-styles]');
    expect(style?.textContent).toContain('.msp-button-button');
    expect(el.classList.contains('msp-root')).toBe(true);

    defineStyles('test/late.css', '.late { color: red }');
    await tick();
    expect(style?.textContent).toContain('.late { color: red }');

    // A nested element is not a root of its own.
    const inner = document.createElement('msp-test-styled') as StyledElement;
    inner.model = 'inner';
    el.append(inner);
    expect(inner.classList.contains('msp-root')).toBe(false);
    // One style element per document.
    expect(document.head.querySelectorAll('style[data-msp-styles]').length).toBe(1);
  });
});
