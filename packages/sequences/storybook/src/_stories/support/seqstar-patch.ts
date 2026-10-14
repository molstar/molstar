/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `storyPatch`: how a Seq* story applies a patch. A story's patch that names the wrong ref, or sets a param the node does
// not have, is not an error the viewer throws: its result carries issues, and a story that drops them looks dead (three
// stories did). This applies the patch and keeps what came back (warnings and errors) in a per-context list that
// `PatchBanner` shows under the viewer, dismissible, and logs it with `console.warn`.

import type { Patch } from "@molstar/seq-core/format";
import type { Issue } from "@molstar/seq-core/model";
import type { Origin, SeqViewerContext } from "@molstar/seq-viewer";

/** What one patch gave back: its label and the issues of severity warning or error. */
export interface PatchNotice {
  readonly id: number;
  readonly label: string;
  readonly issues: readonly Issue[];
}

/** A context's notices: a list the banner reads, changed by `report` and `dismiss`. */
export interface PatchNotices {
  readonly list: readonly PatchNotice[];
  subscribe(fn: () => void): () => void;
  report(label: string, issues: readonly Issue[]): void;
  dismiss(id: number): void;
  clear(): void;
}

const stores = new WeakMap<object, PatchNotices>();
/** The notices of a context (made on first use). */
export function patchNotices(ctx: object): PatchNotices {
  let s = stores.get(ctx);
  if (s) return s;
  let list: readonly PatchNotice[] = [];
  let next = 1;
  const subs = new Set<() => void>();
  const set = (l: readonly PatchNotice[]) => {
    list = l;
    for (const f of subs) f();
  };
  const store: PatchNotices = {
    get list() {
      return list;
    },
    subscribe(fn) {
      subs.add(fn);
      return () => void subs.delete(fn);
    },
    report(label, issues) {
      const bad = issues.filter((i) => i.severity === "warning" || i.severity === "error");
      if (bad.length === 0) return;
      console.warn(`[story patch] "${label}": ${bad.map((i) => `${i.code} ${i.path}: ${i.message}`).join(" | ")}`);
      set([...list, { id: next++, label, issues: bad }]);
    },
    dismiss: (id) => set(list.filter((n) => n.id !== id)),
    clear: () => set([]),
  };
  stores.set(ctx, store);
  s = store;
  return s;
}

/** `ctx.applyPatch`, its issues (warning or error) kept for the banner and logged; the result is as `applyPatch`'s. */
export async function storyPatch(
  ctx: SeqViewerContext,
  patch: Patch,
  label = "patch",
  opts?: { layer?: "base" | "local" | "preview"; snapshot?: string | null; origin?: Origin },
) {
  const r = await ctx.applyPatch(patch, opts);
  patchNotices(ctx).report(label, r.issues);
  return r;
}

/** `ctx.document.propose` (a patch routed by the document port, one undo step), its issues kept and logged. */
export async function storyPropose(ctx: SeqViewerContext, patch: Patch, label: string, origin: Origin) {
  const r = await ctx.document.propose(patch, label, origin);
  patchNotices(ctx).report(label, r.issues);
  return r;
}
