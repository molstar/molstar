/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The remark/rehype adapter (seqspec:libraries#links): for a host whose renderer runs unified (as streamdown does). No
// markdown parsing here; it works on the hast tree the host already has, with no dependency of its own.
//   - `rehypeSeqLinks` (before the sanitizer) marks each `!` link with `data-seq-link`, its id: the occurrence of its
//     href in the document and the href (`seq-link-0:!seq-view=v1`), so a re-rendered, changed document keeps the id of
//     every link it kept, and a link hovered across the change still leaves its own channel. It writes any `:` of the
//     target as `%3A` (a sanitizer reads `!a=b:c` as the protocol `!a=b`; the codec decodes it back). That cannot let
//     a scheme through: a URL whose first character is `!` has none (a scheme starts with a letter), so the link stays
//     relative. Every other link is left as it is, so `javascript:` and `data:` links are still the sanitizer's to drop.
//   - `seqLinkSchema` lets the sanitizer's schema keep the mark.
// The host's link component calls the runner with `data-seq-link` as the link id and its `href`, and `release`s an id
// that leaves the rendered tree.

interface HastNode {
  type: string;
  tagName?: string;
  properties?: Record<string, unknown>;
  children?: HastNode[];
}

export function rehypeSeqLinks(opts: { prefix?: string } = {}) {
  const prefix = opts.prefix ?? "seq-link-";
  return (tree: HastNode) => {
    const seen = new Map<string, number>();
    const visit = (node: HastNode) => {
      const href = node.properties?.href;
      if (node.type === "element" && node.tagName === "a" && typeof href === "string" && href.startsWith("!")) {
        const safe = href.replaceAll(":", "%3A");
        const k = seen.get(safe) ?? 0;
        seen.set(safe, k + 1);
        node.properties = { ...node.properties, href: safe, dataSeqLink: `${prefix}${k}:${safe}` };
      }
      for (const c of node.children ?? []) visit(c);
    };
    visit(tree);
  };
}

/** A hast-util-sanitize schema that also keeps `data-seq-link` on links: `seqLinkSchema(defaultSchema)`. */
export function seqLinkSchema<S extends { attributes?: Record<string, unknown[]> | null }>(base: S): S {
  const a = base.attributes ?? {};
  return { ...base, attributes: { ...a, a: [...(a.a ?? []), "dataSeqLink"] } };
}
