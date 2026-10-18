/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { classes } from '../../components/button/_generated/button.module.css.js';
import { defineElements, MspElement } from '../../dom/element.js';
import { adoptStyles, defineStyles, releaseStyles, styleText } from '../registry.js';

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

describe('styles in shadow roots', () => {
  // A constructed sheet where the environment has them, else the <style> fallback.
  const cssOf = (root: ShadowRoot) => {
    const sheets = (root as ShadowRoot & { adoptedStyleSheets?: CSSStyleSheet[] }).adoptedStyleSheets ?? [];
    const adopted = sheets.map((s) => [...s.cssRules].map((r) => r.cssText).join('\n'));
    const styles = [...root.querySelectorAll('style[data-msp-styles]')].map((s) => s.textContent ?? '');
    return [...adopted, ...styles].join('\n');
  };
  const shadow = () => {
    const host = document.createElement('div');
    document.body.append(host);
    return host.attachShadow({ mode: 'open' });
  };

  it('adopts into several roots, reaches both with a later registration, and releases one', async () => {
    const a = shadow();
    const b = shadow();
    adoptStyles(a);
    adoptStyles(b);
    await tick();
    expect(cssOf(a)).toContain('.msp-button-button');
    expect(cssOf(b)).toContain('.msp-button-button');

    defineStyles('test/shadow-late.css', '.shadow-late { color: red }');
    await tick();
    expect(cssOf(a)).toContain('shadow-late');
    expect(cssOf(b)).toContain('shadow-late');

    releaseStyles(a);
    expect(cssOf(a)).toBe('');
    expect(cssOf(b)).toContain('shadow-late');

    defineStyles('test/shadow-later.css', '.shadow-later { color: blue }');
    await tick();
    expect(cssOf(a)).toBe('');
    expect(cssOf(b)).toContain('shadow-later');

    // Adopting the same root again after the release works.
    adoptStyles(a);
    await tick();
    expect(cssOf(a)).toContain('shadow-later');
    expect(cssOf(a)).toContain('.msp-button-button');
    releaseStyles(a);
    releaseStyles(b);
    expect(cssOf(b)).toBe('');
  });

  it('a release with an update queued does not bring the styles back', async () => {
    const a = shadow();
    adoptStyles(a);
    defineStyles('test/shadow-queued.css', '.shadow-queued { color: red }');
    releaseStyles(a);
    await tick();
    expect(cssOf(a)).toBe('');
  });

  it('keeps the styles of a root adopted twice until it is released twice', async () => {
    const a = shadow();
    adoptStyles(a);
    adoptStyles(a);
    await tick();
    releaseStyles(a);
    expect(cssOf(a)).toContain('.msp-button-button');
    defineStyles('test/shadow-shared.css', '.shadow-shared { color: red }');
    await tick();
    expect(cssOf(a)).toContain('shadow-shared');
    releaseStyles(a);
    expect(cssOf(a)).toBe('');
  });

  it('shares one constructed sheet between the roots of a document, parsed once per registration', async () => {
    // jsdom has constructed sheets but no `adoptedStyleSheets`: give these roots one.
    const withSheets = () => {
      const root = shadow();
      Object.defineProperty(root, 'adoptedStyleSheets', { value: [], writable: true });
      return root as ShadowRoot & { adoptedStyleSheets: CSSStyleSheet[] };
    };
    const a = withSheets();
    const b = withSheets();
    adoptStyles(a);
    adoptStyles(b);
    expect(a.adoptedStyleSheets.length).toBe(1);
    expect(b.adoptedStyleSheets[0]).toBe(a.adoptedStyleSheets[0]);

    const replace = vi.spyOn(CSSStyleSheet.prototype, 'replaceSync');
    defineStyles('test/shared-late.css', '.shared-late { color: red }');
    await tick();
    expect(replace).toHaveBeenCalledTimes(1);
    expect(cssOf(b)).toContain('shared-late');
    replace.mockRestore();

    releaseStyles(a);
    expect(a.adoptedStyleSheets.length).toBe(0);
    expect(cssOf(b)).toContain('shared-late');
    releaseStyles(b);
    expect(b.adoptedStyleSheets.length).toBe(0);
  });

  it('an element releases the styles of the shadow root it leaves, and keeps them when moved within it', async () => {
    const root = shadow();
    const el = document.createElement('msp-test-styled') as StyledElement;
    el.model = 'in a shadow root';
    root.append(el);
    expect(cssOf(root)).toContain('.msp-button-button');

    const box = document.createElement('div');
    root.append(box);
    box.append(el);
    await tick();
    expect(cssOf(root)).toContain('.msp-button-button');

    el.remove();
    await tick();
    expect(cssOf(root)).toBe('');
  });
});
