/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Markdown to a sanitized hast tree with the command links marked: what a host renderer (a Storybook story, Stories'
// streamdown) does with the links adapter, in the form a story can import (the packages are this package's, so a
// story reaches them through this file). Not part of the package.

import rehypeSanitize, { defaultSchema } from "rehype-sanitize";
import remarkParse from "remark-parse";
import remarkRehype from "remark-rehype";
import { unified } from "unified";
import { rehypeSeqLinks, seqLinkSchema } from "../../src/links/remark.ts";

export interface HastNode {
  type: string;
  tagName?: string;
  value?: string;
  properties?: Record<string, unknown>;
  children?: HastNode[];
}

/**
 * `markdown` parsed, `!` links marked with `data-seq-link` (ids from `prefix`, so the same link in two documents lights two
 * channels) and the tree sanitized with the schema that keeps the mark.
 */
export async function markdownToHast(markdown: string, prefix = "seq-link-"): Promise<HastNode> {
  const pipeline = unified()
    .use(remarkParse)
    .use(remarkRehype)
    .use(rehypeSeqLinks, { prefix })
    .use(rehypeSanitize, seqLinkSchema(defaultSchema));
  return (await pipeline.run(pipeline.parse(markdown))) as HastNode;
}
