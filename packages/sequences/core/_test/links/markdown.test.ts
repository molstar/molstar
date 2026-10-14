/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The story helper (markdown.ts): headings and emphasis survive the sanitizer, `!` links are marked with their ids and the
// prefix, other links stay as they are and a script link is dropped.

import { describe, expect, it } from 'vitest';
import { formatLink } from '../../src/links/index.js';
import { type HastNode, markdownToHast } from './markdown.js';

const find = (n: HastNode, tag: string): HastNode[] => [
  ...(n.tagName === tag ? [n] : []),
  ...(n.children ?? []).flatMap((c) => find(c, tag)),
];

describe('markdownToHast', () => {
  it("keeps headings and emphasis, marks command links with the document's prefix, drops a script link", async () => {
    const hl = formatLink('seq-highlight', '1iep/A', { beg: '315', end: '315', numbering: 'auth' });
    const tree = await markdownToHast(
      `# Title\n\nA **bold** [T315](${hl}), [site](https://x.test), [bad](javascript:alert(1)).`,
      'd1-',
    );
    expect(find(tree, 'h1').length).toBe(1);
    expect(find(tree, 'strong').length).toBe(1);
    const [cmd, site, bad] = find(tree, 'a').map((a) => a.properties ?? {});
    expect(cmd.dataSeqLink).toBe(`d1-0:${hl.replaceAll(':', '%3A')}`);
    expect(site).toEqual({ href: 'https://x.test' });
    expect(bad.href).toBeUndefined();
  });
});
