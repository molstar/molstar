/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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
    jest.useFakeTimers();
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
    jest.advanceTimersByTime(600);
    expect(root.hasAttribute('data-scrolling')).toBe(false);
    jest.useRealTimers();
  });
});
