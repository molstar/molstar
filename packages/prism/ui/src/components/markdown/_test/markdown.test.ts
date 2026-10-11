/** @jest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Markdown } from '../markdown.js';

describe('Markdown', () => {
  it('renders headings, paragraphs, lists and inline marks as DOM', () => {
    const el = Markdown({ source: '### Site\nLine one\nline two\n\n- **bold** and *it*\n- `code`\n\n1. first' });
    expect(el.querySelector('h5')?.textContent).toBe('Site');
    expect(el.querySelector('p')?.textContent).toBe('Line one line two');
    expect([...el.querySelectorAll('ul li')].map((li) => li.innerHTML)).toEqual([
      '<strong>bold</strong> and <em>it</em>',
      '<code>code</code>',
    ]);
    expect(el.querySelector('ol li')?.textContent).toBe('first');
  });

  it('never parses HTML, and routes non-web links to onLink', () => {
    const links: string[] = [];
    const el = Markdown({
      source: '<img src=x onerror=alert(1)> [next](#key) [cmd](command:focus) [web](https://molstar.org)',
      onLink: (href) => links.push(href),
    });
    expect(el.querySelector('img')).toBeNull();
    expect(el.textContent).toContain('<img src=x onerror=alert(1)>');
    const [next, cmd, web] = [...el.querySelectorAll('a')];
    next.click();
    cmd.click();
    expect(links).toEqual(['#key', 'command:focus']);
    expect(cmd.getAttribute('href')).toBe('#');
    expect(web.target).toBe('_blank');
    expect(web.rel).toBe('noopener noreferrer');
  });
});
