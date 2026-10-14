/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Seq*/Query and search/Search › Location sets and bookmarks story: p53's Tracks
// view whose section shows three location sets through `highlights`, a SeqQL query node among them, and keeps
// bookmarks. The panel edits the document only through patches (each one undo step): an entry's style, an entry
// removed or added, a bookmark added from the viewer's selection or removed; stepping is the runtime `bookmark.step`.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import type { Loci } from "@molstar/seq-core/model";
import { bookmarkAt, type SeqViewerContext } from "@molstar/seq-viewer";
import { SeqViewerLayout, SeqViewerProvider, useSeqState } from "@molstar/seq-viewer/react";
import { Fragment, type ReactNode, useMemo, useState } from "react";
import { PrismPage } from "../support/seqstar-prism";
import { TextAction, type TextTab, TextTabs } from "../ui";
import { BOOKMARKS, locationSetsDoc, LOCATION_SECTION as SECTION } from "./seqstar-docs-views";
import f from "./seqstar-format.module.css";
import { storyPropose } from "./seqstar-patch";
import { PatchBanner } from "./seqstar-patch-banner";
import { PrismCanvas } from "./seqstar-prism";
import q from "./seqstar-search.module.css";
import { useCaseContext, useStoryUi } from "./seqstar-viewer";

type Json = { kind: string; ref?: string; params?: Record<string, unknown>; children?: Json[] };
type Entry = { data: string; items?: string[] | null; style?: "band" | "highlight" | "accent" };
type Row = { target: string; start: number; stop: number; id: string; label?: string };

/** What the panel offers to add as a location set: the case's annotations and a query node. */
const DATA: readonly { ref: string; label: string }[] = [
  { ref: "cysteines", label: "Cysteines (a SeqQL query)" },
  { ref: "zinc", label: "Zinc ligands" },
  { ref: "dna-contacts", label: "DNA-binding residues" },
  { ref: "domains", label: "Domains" },
  { ref: BOOKMARKS, label: "Bookmarks" },
];
const STYLES: TextTab[] = [
  { value: "band", label: "band", hint: "A full-height band behind the tracks, drawn in the main pass." },
  { value: "highlight", label: "highlight", hint: "The selection's marks on the overlay, in the highlight color." },
  { value: "accent", label: "accent", hint: "The overlay's accent: the highlight color at a quarter of its strength." },
];

const nodeOf = (doc: SeqViewSpecDoc | null, ref: string): Json | undefined => {
  const walk = (n: Json): Json | undefined =>
    n.ref === ref ? n : (n.children ?? []).map(walk).find((x) => x !== undefined);
  return doc ? walk(doc.root as Json) : undefined;
};
const propose = (ctx: SeqViewerContext, title: string, highlights: Entry[]) =>
  void storyPropose(
    ctx,
    { format: "seq_view_spec_patch", ops: [{ op: "set", ref: SECTION, params: { highlights } }] },
    title,
    "host:story",
  );
/** The selection's first part a bookmark can keep, or null. */
function selected(ctx: SeqViewerContext): Loci | null {
  const sel = ctx.selection.current.value;
  return (sel.seq ?? sel.items ?? sel.rows ?? sel.columns ?? null) as Loci | null;
}

/** The section's `highlights` as written, each entry's style switchable and removable, and a picker to add one. */
function SetsPanel({ ctx, entries }: { ctx: SeqViewerContext; entries: Entry[] }) {
  const shown = new Set(entries.map((e) => e.data));
  const restyle = (k: number, style: string) =>
    propose(
      ctx,
      "Style",
      entries.map((e, i) => (i === k ? { ...e, style: style as Entry["style"] } : e)),
    );
  return (
    <div className={q.panel}>
      <span className={q.caption}>
        <code>section.highlights</code> — a channel each, <code>set:{SECTION}/n</code>
      </span>
      {entries.map((e, k) => (
        <div className={q.box} key={e.data}>
          <code className={q.count}>{e.data}</code>
          <TextTabs tabs={STYLES} value={e.style ?? "band"} onValueChange={(v) => restyle(k, v)} />
          <TextAction
            tone="quiet"
            onClick={() =>
              propose(
                ctx,
                "Remove",
                entries.filter((_, i) => i !== k),
              )
            }
          >
            Remove
          </TextAction>
        </div>
      ))}
      <div className={q.box}>
        <span className={q.caption}>Add</span>
        {DATA.filter((d) => !shown.has(d.ref)).map((d) => (
          <TextAction key={d.ref} tone="quiet" onClick={() => propose(ctx, "Add", [...entries, { data: d.ref }])}>
            {d.label}
          </TextAction>
        ))}
      </div>
    </div>
  );
}

/** The bookmarks: the set's items as chips (the current one marked), stepping, adding the selection, removing one. */
function BookmarksPanel({ ctx, rows }: { ctx: SeqViewerContext; rows: Row[] }) {
  const [label, setLabel] = useState("");
  // Re-rendered on every channel change: the step's accent and the selection's light both pass through here.
  useSeqState(ctx.highlight.channels);
  useSeqState(ctx.selection.current);
  const loci = selected(ctx);
  const target = loci && { section: SECTION, loci, label: label || undefined };
  const at = bookmarkAt(ctx, SECTION);
  const ids = [...new Set(rows.map((r) => r.id))];
  const step = (dir: 1 | -1) => ctx.actions.run("bookmark.step", { section: SECTION, dir }, "host:story");
  return (
    <div className={q.panel}>
      <span className={q.caption}>
        Bookmarks — the inline ranges annotation <code>{BOOKMARKS}</code>
      </span>
      <div className={q.box}>
        <span className={q.caption}>Step</span>
        <TextAction onClick={() => step(-1)} disabled={ids.length === 0}>
          Previous
        </TextAction>
        <TextAction onClick={() => step(1)} disabled={ids.length === 0}>
          Next
        </TextAction>
        <span className={q.count} aria-live="polite">
          {at && ids.includes(at) ? `${ids.indexOf(at) + 1} of ${ids.length}` : `${ids.length} bookmarks`}
        </span>
      </div>
      <div className={q.box}>
        <span className={q.caption}>Bookmark the selection as</span>
        <input
          className={q.input}
          value={label}
          placeholder="label (optional)"
          aria-label="Bookmark label"
          onChange={(e) => setLabel(e.target.value)}
        />
        <TextAction
          onClick={() => {
            if (target && ctx.actions.run("bookmark.add", target, "host:story")) setLabel("");
          }}
          title="Select a stretch in the viewer first (drag over the letters)"
          disabled={!target || !ctx.actions.can("bookmark.add", target)}
        >
          Add
        </TextAction>
        <span className={q.caption}>Select a stretch in the viewer first (drag over the letters).</span>
      </div>
      <div className={q.hitList}>
        {ids.map((id) => {
          const own = rows.filter((r) => r.id === id);
          return (
            <Fragment key={id}>
              <button
                type="button"
                className={q.hit}
                aria-current={id === at}
                onClick={() => ctx.actions.run("bookmark.go", { section: SECTION, item: id }, "host:story")}
                // The current one again, or Escape, unselects it: no accent, the focus stays.
                onKeyDown={(e) => e.key === "Escape" && ctx.results.clearCurrent(`bookmarks:${SECTION}`)}
              >
                <span className={q.at}>{id}</span>
                <span className={q.letters}>{own[0].label ?? ""}</span>
                <span className={q.at}>{own.map((r) => `${r.start + 1}–${r.stop}`).join(", ")}</span>
              </button>
              <TextAction
                tone="quiet"
                title={`Remove ${id}`}
                onClick={() => ctx.actions.run("bookmark.remove", { section: SECTION, item: id }, "host:story")}
              >
                ×
              </TextAction>
            </Fragment>
          );
        })}
      </div>
    </div>
  );
}

function Panels({ ctx }: { ctx: SeqViewerContext }) {
  const doc = useSeqState(ctx.document.shown);
  const entries = (nodeOf(doc, SECTION)?.params?.highlights ?? []) as Entry[];
  const rows = (nodeOf(doc, BOOKMARKS)?.params?.value ?? []) as Row[];
  useSeqState(ctx.history.canUndo);
  return (
    <div className={q.controls}>
      <SetsPanel ctx={ctx} entries={entries} />
      <BookmarksPanel ctx={ctx} rows={rows} />
      <TextAction
        tone="quiet"
        onClick={() => ctx.actions.run("history.undo", null)}
        disabled={!ctx.history.canUndo.value}
      >
        Undo
      </TextAction>
    </div>
  );
}

/** The story: the viewer on p53's Tracks view, the two panels above it. */
export function LocationSetsStory(): ReactNode {
  const doc = useMemo(locationSetsDoc, []);
  const m = useCaseContext(doc, { view: "tracks" });
  const ui = useStoryUi(m);
  return (
    <PrismPage
      title="Location sets and bookmarks"
      caption={
        <>
          A location set is a named list of locations: a <code>ranges</code> annotation, or a query node, which is one.
          A section shows its sets persistently through <code>highlights</code>, each entry its own highlight channel (
          <code>band</code> behind the tracks, or the overlay's <code>highlight</code> and <code>accent</code>).
          Bookmarks are such a set the viewer's actions edit: select a stretch in the viewer and bookmark it, step with
          Previous and Next (the current one focused and accented; its pill again, or Escape, unselects it), remove one.
          Every change is a patch, one undo step.
        </>
      }
    >
      <section className={f.card}>
        {m ? (
          <SeqViewerProvider context={m.ctx} ui={ui}>
            <Panels ctx={m.ctx} />
            <PrismCanvas>
              <SeqViewerLayout />
            </PrismCanvas>
            <PatchBanner />
          </SeqViewerProvider>
        ) : (
          <div className={f.status}>Loading…</div>
        )}
      </section>
    </PrismPage>
  );
}
