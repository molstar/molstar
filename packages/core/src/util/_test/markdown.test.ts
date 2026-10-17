/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  isSafeMarkdownHref,
  type MarkdownInline,
  markdownInlineText,
  parseMarkdown,
  parseMarkdownInline,
} from '../markdown.js';

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

  it('reads a link scheme as a browser does: controls and spaces around it stripped, tabs and newlines removed', () => {
    expect(isSafeMarkdownHref('\x01javascript:alert(1)')).toBe(false);
    expect(isSafeMarkdownHref('\x00\x1f javascript:alert(1)\x7f')).toBe(false);
    expect(isSafeMarkdownHref('java\tscript:alert(1)')).toBe(false);
    expect(isSafeMarkdownHref('java\nscript:alert(1)')).toBe(false);
    expect(isSafeMarkdownHref('\x01https://molstar.org')).toBe(true);
    // In the tree: the whole destination is the href, and it is not safe.
    const [link] = parseMarkdownInline('[x](\x01javascript:alert(1))');
    expect(link).toEqual({ kind: 'link', href: '\x01javascript:alert(1)', children: [{ kind: 'text', text: 'x' }] });
    expect(isSafeMarkdownHref((link as { href: string }).href)).toBe(false);
    // An entity-encoded scheme is decoded first, so the check sees it.
    const [encoded] = parseMarkdownInline('[x](&#106;avascript:alert(1))') as { href: string }[];
    expect(encoded.href).toBe('javascript:alert(1)');
    expect(isSafeMarkdownHref(encoded.href)).toBe(false);
  });
});

describe('links', () => {
  it('keeps balanced parentheses in a destination and decodes its entities and escapes', () => {
    expect(parseMarkdownInline('[w](https://en.wikipedia.org/wiki/Foo_(bar)) after')).toEqual([
      { kind: 'link', href: 'https://en.wikipedia.org/wiki/Foo_(bar)', children: [{ kind: 'text', text: 'w' }] },
      { kind: 'text', text: ' after' },
    ]);
    expect(parseMarkdownInline('[a](x?b=1&amp;c=2 "t")')).toEqual([
      { kind: 'link', href: 'x?b=1&c=2', children: [{ kind: 'text', text: 'a' }] },
    ]);
    expect(parseMarkdownInline('[a](x\\)y)')).toEqual([
      { kind: 'link', href: 'x)y', children: [{ kind: 'text', text: 'a' }] },
    ]);
    // An unclosed parenthesis is no link.
    expect(markdownInlineText(parseMarkdownInline('[a](x(y)'))).toBe('[a](x(y)');
  });
});

describe('entities', () => {
  it('decodes only the known names, never what an object inherits', () => {
    expect(parseMarkdownInline('&constructor; &toString; &valueOf; &hasOwnProperty;')).toEqual([
      { kind: 'text', text: '&constructor; &toString; &valueOf; &hasOwnProperty;' },
    ]);
    expect(parseMarkdownInline('&lt;&nbsp;&#65;&#x42;')).toEqual([{ kind: 'text', text: '<\u00a0AB' }]);
  });
});

describe('nesting', () => {
  it('caps the depth: deeper marks stay plain text, and nothing overflows the stack', () => {
    const deep = `${'**'.repeat(30000)}a${'**'.repeat(30000)}`;
    let nodes: readonly MarkdownInline[] = [];
    expect(() => (nodes = parseMarkdownInline(deep))).not.toThrow();
    let depth = 0;
    for (let n = nodes[0]; n && n.kind === 'strong'; n = n.children[0]) depth++;
    expect(depth).toBeGreaterThan(1);
    expect(depth).toBeLessThanOrEqual(8);
    expect(() => parseMarkdownInline(`${'['.repeat(30000)}a${'](x)'.repeat(30000)}`)).not.toThrow();
    expect(() => parseMarkdownInline(`${'*'.repeat(60000)}a${'*'.repeat(60000)}`)).not.toThrow();
    // Ordinary nesting is still read.
    expect(parseMarkdownInline('[**a *b* c**](x)')).toEqual([
      {
        kind: 'link',
        href: 'x',
        children: [
          {
            kind: 'strong',
            children: [
              { kind: 'text', text: 'a ' },
              { kind: 'em', children: [{ kind: 'text', text: 'b' }] },
              { kind: 'text', text: ' c' },
            ],
          },
        ],
      },
    ]);
  });
});

describe('adversarial input', () => {
  // Each was quadratic or worse (a heading with a run of spaces took 72 s at 8 KB); linear, 100 KB parses in a few ms.
  // The bound is generous (an order of magnitude over a slow machine), far under what a quadratic needs.
  const size = 100_000;
  const fill = (unit: string) => unit.repeat(Math.ceil(size / unit.length)).slice(0, size);
  const inputs: Record<string, string> = {
    'heading with a run of spaces': `# a${' '.repeat(size)}b`,
    'heading with a run of hashes': `# ${'#'.repeat(size)}a`,
    'list item with a line separator': `- ${' '.repeat(size)}a\u2028`,
    'unclosed strong': fill('**a '),
    'unclosed strong underscores': fill('__a '),
    'unclosed emphasis': fill('*a '),
    'unclosed underscore emphasis': fill('_a '),
    'unclosed brackets': fill('[a '),
    'unclosed destinations': fill('[a]('),
    'unclosed parentheses': fill('[a](('),
    'unclosed titles': fill('[a](b "'),
    'unclosed images': fill('![a]('),
    'unclosed tags': fill('<a '),
    'unclosed autolinks': fill('<http:'),
    'unclosed code': fill('`a ``b '),
    'a run of spaces': `a${' '.repeat(size)}b`,
    'deep strong': `${'**'.repeat(size / 4)}a${'**'.repeat(size / 4)}`,
    'mixed openers': fill('*_[`<!**__'),
  };
  for (const [name, text] of Object.entries(inputs)) {
    it(name, () => {
      const start = performance.now();
      parseMarkdown(text);
      expect(performance.now() - start).toBeLessThan(200);
    });
  }
});
