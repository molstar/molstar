/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import rehypeSanitize, { defaultSchema } from "rehype-sanitize";
import remarkParse from "remark-parse";
import remarkRehype from "remark-rehype";
import { unified } from "unified";
import { describe, expect, it } from "vitest";
import { createLinkRunner, formatLink, parseLink } from "../../src/links/index.js";
import { rehypeSeqLinks, seqLinkSchema } from "../../src/links/remark.js";

interface Node {
  type: string;
  tagName?: string;
  properties?: Record<string, unknown>;
  children?: Node[];
}
/** The `a` elements of a hast tree, in order. */
const links = (n: Node): Node[] => [...(n.tagName === "a" ? [n] : []), ...(n.children ?? []).flatMap(links)];

async function render(md: string) {
  const tree = unified()
    .use(remarkParse)
    .use(remarkRehype)
    .use(rehypeSeqLinks)
    .use(rehypeSanitize, seqLinkSchema(defaultSchema));
  return links((await tree.run(tree.parse(md))) as Node).map((a) => a.properties ?? {});
}

describe("remark adapter", () => {
  const sel = formatLink("seq-select", "unit:1,unit:2", { start: 3, stop: 9 });
  const md = `[a](${sel}) [b](https://example.com/x) [c](javascript:alert(1)) [d](${formatLink("seq-view", "v1")})`;

  it("marks ! links, leaves https links alone, and the sanitizer drops javascript:", async () => {
    const [a, b, c, d] = await render(md);
    expect(a.dataSeqLink).toBe(`seq-link-0:${sel}`);
    expect(d.dataSeqLink).toBe("seq-link-0:!seq-view=v1");
    expect(b).toEqual({ href: "https://example.com/x" });
    expect(c.href).toBeUndefined();
    expect(c.dataSeqLink).toBeUndefined();
  });
  it("keeps the sanitizer's scheme rules: data:, javascript: in any case, and a ! link with a scheme inside", async () => {
    const out = await render(
      "[e](data:text/html;base64,PHNjcmlwdD4=) [f](JavaScript:alert(1)) [g](vbscript:x) [h](!javascript:alert(1))",
    );
    expect(out.slice(0, 3).map((p) => p.href)).toEqual([undefined, undefined, undefined]);
    // A ! link is relative (no scheme starts with `!`) and has no colon left; its key is nobody's.
    expect(out[3].href).toBe("!javascript%3Aalert(1)");
    expect(createLinkRunner(() => {}).run(out[3].href as string, "click", "h")).toBe(false);
  });
  it("keeps a target with a colon through the sanitizer, formatted or written by hand, and it still runs", async () => {
    const [a] = await render(md);
    expect(a.href).toBeTypeOf("string");
    const calls: unknown[] = [];
    const runner = createLinkRunner((inv) => calls.push(inv));
    expect(runner.run(a.href as string, "click", a.dataSeqLink as string)).toBe(true);
    const [e] = (await render("[e](!seq-select=unit:3)")) as { href: string }[];
    expect(runner.run(e.href, "click", "e")).toBe(true);
    expect(calls).toEqual([
      {
        command: "seqv:select",
        params: {
          target: { refs: ["unit:1", "unit:2"], args: { start: 3, stop: 9 } },
          mode: "set",
          closure: "containing",
        },
      },
      {
        command: "seqv:select",
        params: { target: { refs: ["unit:3"], args: {} }, mode: "set", closure: "containing" },
      },
    ]);
  });
  it("a formatted ref with parentheses stays one link", async () => {
    const href = formatLink("seq-view", "v(1))");
    const [a] = await render(`[a](${href}) after`);
    expect(parseLink(a.href as string)[0].value).toBe("v(1))");
  });
  it("a sanitizer without the schema helper drops the mark, which is why the helper exists", async () => {
    const tree = unified().use(remarkParse).use(remarkRehype).use(rehypeSeqLinks).use(rehypeSanitize);
    const [a] = links((await tree.run(tree.parse(md))) as Node);
    expect(a.properties?.dataSeqLink).toBeUndefined();
  });
  it("ids are the href and its occurrence, kept when the document changes around a link", async () => {
    const tree = unified().use(remarkParse).use(remarkRehype).use(rehypeSeqLinks, { prefix: "p-" });
    const ids = async (doc: string) =>
      links((await tree.run(tree.parse(doc))) as Node).map((a) => a.properties?.dataSeqLink as string);
    expect(await ids("[a](!seq-view=a) [b](!seq-view=b) [b again](!seq-view=b)")).toEqual([
      "p-0:!seq-view=a",
      "p-0:!seq-view=b",
      "p-1:!seq-view=b",
    ]);
    // Hovered b, then a link is written before it: its leave still clears the channel its enter lit.
    const r: { channel: string; target: unknown }[] = [];
    const runner = createLinkRunner((inv) => r.push(inv.params as (typeof r)[number]));
    const [, b] = await ids("[a](!seq-highlight=a) [b](!seq-highlight=b)");
    runner.run("!seq-highlight=b", "enter", b);
    const [, , b2] = await ids("[x](!seq-highlight=x) [a](!seq-highlight=a) [b](!seq-highlight=b)");
    runner.run("!seq-highlight=b", "leave", b2);
    expect(r.map((x) => [x.channel, x.target === null])).toEqual([
      ["links:p-0:!seq-highlight=b", false],
      ["links:p-0:!seq-highlight=b", true],
    ]);
  });
});
