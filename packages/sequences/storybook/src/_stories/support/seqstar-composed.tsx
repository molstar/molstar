/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Seq*/Cases/Default UIs › composed stories: the viewer's chrome as components a host places
// (seqspec:interaction#slots-and-chrome). The viewer is `clean` (the canvas, hover, the tooltip); the host puts the view
// tabs, each section's switches, the export menu and the issue badge in a header row of its own above it, and the legend
// below it, all against the same context. React through the provider; the DOM through the element's `mountedUI`.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { shownSections } from "@molstar/seq-viewer";
import {
  browserExportAction,
  Export as DomExport,
  IssueBadge as DomIssueBadge,
  Legend as DomLegend,
  ViewTabs as DomViewTabs,
  defineSeqViewer,
  downloadSvg,
  type SeqViewerElement,
} from "@molstar/seq-viewer/dom";
import {
  Export,
  IssueBadge,
  Legend,
  SectionSwitches,
  SeqViewerLayout,
  SeqViewerProvider,
  useDerived,
  ViewTabs,
} from "@molstar/seq-viewer/react";
import { type ReactNode, useEffect, useRef } from "react";
import "@molstar/seq-viewer/skin.css";
import { PrismCanvas } from "./seqstar-prism";
import p from "./seqstar-prism.module.css";
import { type Mounted, useCaseContext, useStoryUi } from "./seqstar-viewer";

/** The browser's SVG export, which `SeqViewer` and the element register themselves; a host's context registers it. */
const EXPORT = [{ actions: [browserExportAction({ offer: downloadSvg })] }];

/** The host's header: the tabs, then each shown section's switches, the export menu and the issues. */
function HostHeader({ m }: { m: Mounted }) {
  const ctx = m.ctx;
  const sections = useDerived([ctx.graph.assembled], () => shownSections(ctx), { deps: [ctx] });
  return (
    <div className={p.hostBar}>
      <ViewTabs context={ctx} />
      <span className={p.grow} />
      {sections.map((s) => (
        <SectionSwitches key={s} context={ctx} section={s} />
      ))}
      <Export context={ctx} />
      <IssueBadge context={ctx} />
    </div>
  );
}

/** The React composition: the host's header, the clean viewer, the legend, on a canvas of `color`. */
export function ComposedReact({ doc, color }: { doc: SeqViewSpecDoc; color?: string }): ReactNode {
  const m = useCaseContext(doc, { chrome: "clean", registry: EXPORT });
  const ui = useStoryUi(m);
  if (!m) return <div className={p.caption}>Loading…</div>;
  return (
    <PrismCanvas color={color}>
      <SeqViewerProvider context={m.ctx} ui={ui}>
        <div className={p.composed}>
          <HostHeader m={m} />
          <SeqViewerLayout />
          <Legend context={m.ctx} />
        </div>
      </SeqViewerProvider>
    </PrismCanvas>
  );
}

/**
 * The DOM composition: a `clean` `<seq-viewer>`, and the tabs, export, issues and legend mounted by the host into its
 * own elements around it, sharing the element's UI (`mountedUI`: its hints, menus and surfaces). Parts outside the
 * element's shadow root take the skin from the page (`@molstar/seq-viewer/skin.css`).
 */
export function ComposedDom({ make }: { make: (el: SeqViewerElement) => void }): ReactNode {
  const top = useRef<HTMLDivElement>(null);
  const middle = useRef<HTMLDivElement>(null);
  const bottom = useRef<HTMLDivElement>(null);
  useEffect(() => {
    defineSeqViewer();
    const el = document.createElement("seq-viewer");
    el.chrome = "clean";
    make(el);
    middle.current?.append(el);
    const ctx = el.context;
    const ui = el.mountedUI ?? undefined;
    const [bar, below] = [top.current, bottom.current];
    if (!ctx || !bar || !below) return () => el.remove();
    const grow = document.createElement("span");
    grow.className = p.grow;
    const parts = [DomViewTabs(bar, ctx, { ui })];
    bar.append(grow);
    parts.push(DomExport(bar, ctx, { ui }), DomIssueBadge(bar, ctx, { ui }), DomLegend(below, ctx, { ui }));
    return () => {
      for (const part of parts) part.dispose();
      grow.remove();
      el.remove();
    };
  }, [make]);
  return (
    <PrismCanvas>
      <div className={p.composed}>
        <div ref={top} className={p.hostBar} />
        <div ref={middle} />
        <div ref={bottom} />
      </div>
    </PrismCanvas>
  );
}
