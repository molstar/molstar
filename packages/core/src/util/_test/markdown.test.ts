/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { isSafeMarkdownHref, markdownInlineText, parseMarkdown, parseMarkdownInline } from '../markdown.js';

describe('parseMarkdownInline', () => {
  it('reads emphasis, code, links and breaks', () => {
    expect(parseMarkdownInline('a **b** *c* `d` [e](https://x.org) f')).toEqual([
      { kind: 'text', text: 'a ' },
      { kind: 'strong', children: [{ kind: 'text', text: 'b' }] },
      { kind: 'text', text: ' ' },
      { kind: 'em', children: [{ kind: 'text', text: 'c' }] },
      { kind: 'text', text: ' ' },
      { kind: 'code', text: 'd' },
      { kind: 'text', text: ' ' },
      { kind: 'link', href: 'https://x.org', children: [{ kind: 'text', text: 'e' }] },
      { kind: 'text', text: ' f' },
    ]);
    expect(parseMarkdownInline('one  \ntwo\nthree')).toEqual([
      { kind: 'text', text: 'one' },
      { kind: 'break' },
      { kind: 'text', text: 'two three' },
    ]);
  });

  it('keeps underscores inside words, escapes, entities and nesting', () => {
    expect(markdownInlineText(parseMarkdownInline('snake_case_name'))).toBe('snake_case_name');
    expect(parseMarkdownInline('\\*not em\\* &amp; &#x41;')).toEqual([{ kind: 'text', text: '*not em* & A' }]);
    expect(parseMarkdownInline('**bold *and em***')).toEqual([
      {
        kind: 'strong',
        children: [
          { kind: 'text', text: 'bold ' },
          { kind: 'em', children: [{ kind: 'text', text: 'and em' }] },
        ],
      },
    ]);
  });

  it('drops HTML and keeps it out of the tree entirely', () => {
    const nodes = parseMarkdownInline('x<img src=y onerror=alert(1)><script>alert(2)</script> z');
    expect(markdownInlineText(nodes)).toBe('xalert(2) z');
    expect(JSON.stringify(nodes)).not.toMatch(/<|onerror/);
  });

  it('reads autolinks and image alt text', () => {
    expect(parseMarkdownInline('<https://molstar.org> ![a cat](c.png)')).toEqual([
      { kind: 'link', href: 'https://molstar.org', children: [{ kind: 'text', text: 'https://molstar.org' }] },
      { kind: 'text', text: ' a cat' },
    ]);
  });
});

describe('parseMarkdown', () => {
  it('reads headings, paragraphs, lists, code and rules', () => {
    const blocks = parseMarkdown(
      '# Title\n\nSome text\nwrapped.\n\n- one\n- two\n  - nested\n1. first\n\n```js\nlet x = 1;\n```\n\n---',
    );
    expect(blocks.map((b) => b.kind)).toEqual(['heading', 'paragraph', 'list', 'list', 'code', 'rule']);
    expect(blocks[0]).toEqual({ kind: 'heading', level: 1, children: [{ kind: 'text', text: 'Title' }] });
    expect(blocks[1]).toEqual({ kind: 'paragraph', children: [{ kind: 'text', text: 'Some text wrapped.' }] });
    const bullets = blocks[2] as Extract<(typeof blocks)[number], { kind: 'list' }>;
    expect(bullets.ordered).toBe(false);
    expect(bullets.items.length).toBe(2);
    expect(bullets.items[1].map((b) => b.kind)).toEqual(['paragraph', 'list']);
    const ordered = blocks[3] as Extract<(typeof blocks)[number], { kind: 'list' }>;
    expect([ordered.ordered, ordered.start]).toEqual([true, 1]);
    expect(blocks[4]).toEqual({ kind: 'code', lang: 'js', text: 'let x = 1;' });
  });

  it('says which links a renderer may follow', () => {
    expect(isSafeMarkdownHref('https://molstar.org')).toBe(true);
    expect(isSafeMarkdownHref('mailto:a@b.c')).toBe(true);
    expect(isSafeMarkdownHref('docs/page.html')).toBe(true);
    expect(isSafeMarkdownHref('javascript:alert(1)')).toBe(false);
    expect(isSafeMarkdownHref(' JavaScript:alert(1)')).toBe(false);
    expect(isSafeMarkdownHref('data:text/html,x')).toBe(false);
  });
});
