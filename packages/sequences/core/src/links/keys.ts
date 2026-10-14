/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Link keys (seqspec:libraries#links): what each key of a link means as an invocation of a registered command, as plain
// data (seqspec:format#descriptions-and-invocations). The ids are the viewer's built-ins (seqspec:viewer#commands); an
// application's keys map to its own.

import type { Invocation } from "@molstar/seq-core/format";
import type { ParsedLink } from "./codec.js";

export interface LinkKey {
  /** "seq-focus", or an application's "<ns>:<name>". */
  readonly key: string;
  /** What the link's value is, for builders and the reference. */
  readonly value: "refs" | "ref" | "none";
  /** Null: nothing to run for this link. */
  readonly enter?: (l: ParsedLink) => Invocation | null;
  readonly leave?: (l: ParsedLink) => Invocation | null;
  readonly click?: (l: ParsedLink) => Invocation | null;
}

/** Link args arrive as strings: `start` and `stop` as integers, `items` as a list, absent args left out. */
function target(l: ParsedLink) {
  const a = l.args;
  const refs = l.value.split(",").filter((r) => r !== "");
  const args: Record<string, unknown> = {};
  for (const k of ["start", "stop"]) {
    const n = Number(a[k]);
    if (a[k] && Number.isInteger(n)) args[k] = n;
  }
  for (const k of ["beg", "end", "numbering"]) if (a[k]) args[k] = a[k];
  if (a.items) args.items = a.items.split(",").filter((i) => i !== "");
  return refs.length === 0 ? null : { refs, args };
}

const channel = (l: ParsedLink) => `links:${l.link}`;
const highlight = (l: ParsedLink): Invocation | null => {
  const t = target(l);
  return t && { command: "seqv:highlight", params: { channel: channel(l), target: t, closure: "containing" } };
};
const unhighlight = (l: ParsedLink): Invocation => ({
  command: "seqv:highlight",
  params: { channel: channel(l), target: null },
});
const hover = { value: "refs", enter: highlight, leave: unhighlight } as const;

/** The five keys of the format. `section` is `seq-focus`'s alone. */
export const DefaultLinkKeys: readonly LinkKey[] = [
  { key: "seq-highlight", ...hover },
  {
    key: "seq-select",
    ...hover,
    click: (l) => {
      const t = target(l);
      return t && { command: "seqv:select", params: { target: t, mode: "set", closure: "containing" } };
    },
  },
  {
    key: "seq-focus",
    ...hover,
    click: (l) => {
      const t = target(l);
      return (
        t && { command: "seqv:focus", params: l.args.section ? { target: t, section: l.args.section } : { target: t } }
      );
    },
  },
  {
    key: "seq-snapshot",
    value: "ref",
    click: (l) => ({ command: "seqv:setSnapshot", params: { ref: l.value || null } }),
  },
  {
    key: "seq-view",
    value: "ref",
    click: (l) => (l.value ? { command: "seqv:setView", params: { ref: l.value } } : null),
  },
];
