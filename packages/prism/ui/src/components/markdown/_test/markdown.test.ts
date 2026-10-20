/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { Markdown, type MarkdownHost } from '../markdown.js';

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

describe('emphasis', () => {
  it('holds links and code: a command link in bold goes to the host', () => {
    const calls: string[] = [];
    const host: MarkdownHost = {
      link: (href) => (href.startsWith('!') ? (e) => calls.push(`${href} ${e}`) : undefined),
    };
    const el = Markdown({
      source: 'a **[heme group](!query%3Dresn%20HEM%26lang%3Dpymol)** and *[web](https://molstar.org) `x`*',
      host,
    });
    const strongLink = el.querySelector('strong a') as HTMLAnchorElement;
    expect(strongLink.textContent).toBe('heme group');
    strongLink.click();
    expect(calls).toEqual(['!query%3Dresn%20HEM%26lang%3Dpymol click']);
    expect(el.querySelector('em a')?.getAttribute('href')).toBe('https://molstar.org');
    expect(el.querySelector('em code')?.textContent).toBe('x');
    expect(el.textContent).not.toContain('](');
  });
});

describe('character references', () => {
  it('are the characters they name in text, emphasis, headings, links, tables and alt text; not in code', () => {
    const el = Markdown({
      source: [
        '## A&nbsp;title &amp; more',
        '[&#x2B05; Back to *Basic workflow*](#basic) and&nbsp;then &#8212; **&lt;b&gt;** *&quot;it&quot;*',
        '',
        '| a&nbsp;b | &#65; |',
        '|---|---|',
        '| &copy; | ![x&amp;y](img.png) |',
        '',
        '`&nbsp;` &unknownname;',
        '```',
        '&amp; stays',
        '```',
      ].join('\n'),
    });
    expect(el.querySelector('h4')?.textContent).toBe('A\u00a0title & more');
    expect(el.querySelector('a')?.textContent).toBe('\u2B05 Back to Basic workflow');
    expect(el.querySelector('a em')?.textContent).toBe('Basic workflow');
    const p = el.querySelector('p')!;
    expect(p.textContent).toContain('and\u00a0then \u2014 <b> "it"');
    // Decoded markup is text, never elements.
    expect(p.querySelector('b')).toBeNull();
    expect([...el.querySelectorAll('th, td')].map((c) => c.textContent)).toEqual(['a\u00a0b', 'A', '\u00a9', '']);
    expect(el.querySelector('img')?.alt).toBe('x&y');
    const last = el.querySelectorAll('p')[1];
    expect(last.querySelector('code')?.textContent).toBe('&nbsp;');
    expect(last.textContent).toContain('&unknownname;');
    expect(el.querySelector('pre code')?.textContent).toBe('&amp; stays');
  });
});

describe('Markdown with a host', () => {
  it('renders GFM tables with their alignment, fenced code and rules', () => {
    const el = Markdown({
      source: '| Residue | Count |\n|:--|--:|\n| **HEM** | 1 |\n| A\\|B | 2 |\n\n```\nx = 1\n```\n\n---',
    });
    expect([...el.querySelectorAll('th')].map((c) => c.textContent)).toEqual(['Residue', 'Count']);
    expect([...el.querySelectorAll('tbody tr')].map((r) => r.textContent)).toEqual(['HEM1', 'A|B2']);
    expect(el.querySelector('td strong')?.textContent).toBe('HEM');
    expect((el.querySelectorAll('td')[1] as HTMLElement).style.textAlign).toBe('right');
    expect(el.querySelector('pre code')?.textContent).toBe('x = 1');
    expect(el.querySelector('hr')).not.toBeNull();
  });

  it('sends the links the host claims to it, on click, pointer enter and leave; the rest as before', () => {
    const events: string[] = [];
    const links: string[] = [];
    const host: MarkdownHost = {
      link: (href) => (href.startsWith('!') ? (event) => events.push(`${href} ${event}`) : undefined),
    };
    const el = Markdown({
      source: '[focus](!focus-refs=a) [next](#key) [web](https://molstar.org)',
      onLink: (href) => links.push(href),
      host,
    });
    const [focus, next, web] = [...el.querySelectorAll('a')];
    focus.dispatchEvent(new Event('pointerenter'));
    focus.click();
    focus.dispatchEvent(new Event('pointerleave'));
    next.click();
    expect(events).toEqual(['!focus-refs=a enter', '!focus-refs=a click', '!focus-refs=a leave']);
    expect(links).toEqual(['#key']);
    expect(focus.getAttribute('href')).toBe('#');
    expect(web.target).toBe('_blank');
  });

  it('resolves images through the host: a URL, a promise of one, an element, or the source as it is', async () => {
    const swatch = document.createElement('span');
    let resolve!: (url: string) => void;
    const host: MarkdownHost = {
      image: (src) =>
        src === 'a.png'
          ? 'data:image/png;base64,AA'
          : src === 'b.png'
            ? new Promise<string>((r) => {
                resolve = r;
              })
            : src.startsWith('!')
              ? swatch
              : undefined,
    };
    const el = Markdown({ source: '![A](a.png) ![B](b.png) ![](!color-swatch=red) ![C](https://x.org/c.png)', host });
    const [a, b, c] = [...el.querySelectorAll('img')];
    expect(a.getAttribute('src')).toBe('data:image/png;base64,AA');
    expect(a.alt).toBe('A');
    expect(b.hidden).toBe(true);
    resolve('data:image/png;base64,BB');
    await Promise.resolve();
    await Promise.resolve();
    expect(b.getAttribute('src')).toBe('data:image/png;base64,BB');
    expect(b.hidden).toBe(false);
    expect(el.contains(swatch)).toBe(true);
    expect(c.getAttribute('src')).toBe('https://x.org/c.png');
  });

  it("shows the host's media above the text", () => {
    const audio = document.createElement('audio');
    const el = Markdown({ source: 'Listen', host: { media: observableState<Node | null>(audio).asReadonly } });
    expect(el.firstElementChild?.contains(audio)).toBe(true);
  });
});

describe('Markdown in a view', () => {
  it('ends the hover of a link removed under the pointer', () => {
    const events: string[] = [];
    const host: MarkdownHost = { link: () => (event) => events.push(event) };
    const target = document.createElement('div');
    const dispose = render(() => Markdown({ source: '[hl](!highlight-refs=a)', host }), target);
    target.querySelector('a')!.dispatchEvent(new Event('pointerenter'));
    dispose();
    expect(events).toEqual(['enter', 'leave']);
  });
});
