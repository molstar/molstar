/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { ScrollArea } from '../scroll-area.js';

describe('ScrollArea', () => {
  it('puts role, aria, ref and children on the content, the scroll on the viewport, and a bar per axis', () => {
    const container = document.createElement('div');
    let content: HTMLElement | undefined;
    let viewport: HTMLElement | undefined;
    render(
      () => (
        <ScrollArea
          class="area"
          contentClass="list"
          role="listbox"
          aria-label="Items"
          ref={(el) => {
            content = el;
          }}
          viewportRef={(el) => {
            viewport = el;
          }}
        >
          <div>one</div>
        </ScrollArea>
      ),
      container,
    );
    const root = container.firstElementChild as HTMLElement;
    expect(root.classList.contains('area')).toBe(true);
    expect(root.dataset.axis).toBe('y');
    expect(content?.getAttribute('role')).toBe('listbox');
    expect(content?.getAttribute('aria-label')).toBe('Items');
    expect(content?.classList.contains('list')).toBe(true);
    expect(content?.textContent).toBe('one');
    expect(viewport?.firstElementChild).toBe(content);
    expect(root.querySelectorAll('[data-orientation]').length).toBe(1);
    expect(root.querySelector('[data-orientation]')?.getAttribute('aria-hidden')).toBe('true');
  });

  it('shows itself as scrolling for a moment after a scroll', () => {
    vi.useFakeTimers();
    const container = document.createElement('div');
    let viewport!: HTMLElement;
    render(
      () => (
        <ScrollArea
          axis="both"
          viewportRef={(el) => {
            viewport = el;
          }}
        >
          x
        </ScrollArea>
      ),
      container,
    );
    const root = container.firstElementChild as HTMLElement;
    expect(root.querySelectorAll('[data-orientation]').length).toBe(2);
    viewport.dispatchEvent(new Event('scroll'));
    expect(root.hasAttribute('data-scrolling')).toBe(true);
    vi.advanceTimersByTime(600);
    expect(root.hasAttribute('data-scrolling')).toBe(false);
    vi.useRealTimers();
  });

  it('is hidden at rest and shown by hover, scrolling, focus inside or a drag: opacity only, no layout', () => {
    // The rules are CSS: the generated module holds them (a jsdom page has no hover to try).
    const css = readFileSync(path.resolve(import.meta.dirname, '../_generated/scroll-area.module.css.ts'), 'utf8');
    const bar = /\.msp-scroll-area-scrollbar \{[^}]*\}/.exec(css)?.[0] ?? '';
    expect(bar).toContain('opacity: 0');
    expect(bar).toContain('transition: opacity');
    const shown = /([^{}]*)\{\s*opacity: 1;?\s*\}/.exec(css)?.[1] ?? '';
    for (const state of [':hover', ':focus-within', '[data-scrolling]', '[data-dragging]'])
      expect(shown).toContain(state);
    // Showing it changes nothing but opacity: no display, width or position in the shown rule.
    expect(/opacity: 1;?\s*\}/.test(css)).toBe(true);
  });

  it('marks the track while its thumb is dragged', () => {
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ScrollArea>x</ScrollArea>, container);
    const root = container.firstElementChild as HTMLElement;
    const track = root.querySelector('[data-orientation]') as HTMLElement;
    // Nothing overflows in jsdom: a press does nothing (no drag, no mark).
    track.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, button: 0 }));
    expect(track.hasAttribute('data-dragging')).toBe(false);
    container.remove();
  });
});
