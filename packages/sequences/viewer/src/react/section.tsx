/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `ViewSurface` (seqspec:interaction#the-default-ui): the current view's sections in a column, each a frame with its
// header (its row holds the hidden tracks' chips) and a mounted surface, plus the chrome placed over the canvas: the
// row controls and the track header's menu. Those follow the viewport through styles set in subscriptions; React
// re-renders them only when the rows in view or the track
// shown change. Keyboard focus on a surface is shown only when it came from Tab (seqspec:interaction#accessibility). As
// src/dom/section.ts.

import { Height } from '@molstar/seq-core/format';
import { type ComponentType, type RefObject, useEffect, useLayoutEffect, useRef, useState } from 'react';
import { later, off, setAttr, watch } from '../dom/core.js';
import { focusIndicator, placeScroll, setStyle, type TrackBox, trackBox } from '../dom/shared.js';
import { shownSections } from '../index.js';
import {
  type PartProps,
  type RowProps,
  type SectionProps,
  type TrackProps,
  useDerived,
  useLanding,
  useScope,
} from './hooks.js';

const toParent = (el: HTMLElement) => el.parentElement ?? el;

export function ViewSurface({ context: ctx }: PartProps) {
  const sections = useDerived([ctx.graph.assembled], () => shownSections(ctx), { deps: [ctx] });
  return (
    <div className="seqv-view">
      {sections.map((ref) => (
        <SectionFrame key={ref} context={ctx} section={ref} />
      ))}
    </div>
  );
}

function SectionFrame({ context: ctx, section }: SectionProps) {
  const { ui, surfaces } = useScope();
  // The section's chrome, as the layout's composition places it (seqspec:interaction#slots-and-chrome).
  const placed = <K extends 'SectionHeader' | 'TrackActions' | 'RowControl'>(k: K) =>
    ui.placed.has(k) ? ui.components[k] : null;
  const [SectionHeader, TrackActions, RowControl] = [
    placed('SectionHeader'),
    placed('TrackActions'),
    placed('RowControl'),
  ];
  const root = useRef<HTMLElement>(null);
  const host = useRef<HTMLDivElement>(null);
  // Focus inside a section that leaves lands on the viewer's root.
  useLanding(root, toParent);
  useEffect(() => {
    const el = host.current as HTMLDivElement;
    const s = ui.surface(ctx, section, el);
    surfaces.set(section, s.element);
    const subs = [
      ctx.layout.contentHeight(section)?.subscribe((height) => {
        setStyle(el, 'height', `${Math.max(40, Math.ceil(height))}px`);
        setAttr(el, 'data-height', ctx.layout.get(section)?.current.heightMode === Height.box ? 'box' : 'content');
      }),
    ];
    const indicator = focusIndicator(el, s.element);
    return () => {
      indicator();
      off(subs);
      if (surfaces.get(section) === s.element) surfaces.delete(section);
      s.dispose();
    };
  }, [ctx, section, ui.surface, surfaces]);
  return (
    <section ref={root} className="seqv-section">
      {SectionHeader && <SectionHeader context={ctx} section={section} />}
      <div className="seqv-stage">
        <div ref={host} className="seqv-surface" />
        {RowControl && <RowControls context={ctx} section={section} Control={RowControl} />}
        {TrackActions && <TrackActionsAt context={ctx} section={section} Actions={TrackActions} />}
      </div>
    </section>
  );
}

/**
 * The `TrackActions` slot at the header `trackBox` decides (the hovered, else the keyboard cursor's; the rules are
 * shared with the DOM UI), kept while the pointer is over it, while its menu is open and, with nothing hovered, while
 * focus is inside it. React renders only when the track changes; its line (a rows track's hovered row) follows each
 * hover by styles.
 */
function TrackActionsAt({ context: ctx, section, Actions }: SectionProps & { Actions: ComponentType<TrackProps> }) {
  const [track, setTrack] = useState<string | null>(null);
  const [tracker, setTracker] = useState<TrackBox | null>(null);
  const box = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const t = trackBox(ctx, section, { box: () => box.current, show: setTrack });
    setTracker(t);
    const subs = watch([ctx.hover.live.version, ctx.ui.keyCursor], () => t.update());
    return () => {
      off(subs);
      t.dispose();
      setTracker(null);
      setTrack(null);
    };
  }, [ctx, section]);
  if (!track || !tracker) return null;
  return (
    <TrackBoxView
      key={track}
      context={ctx}
      section={section}
      track={track}
      box={box}
      tracker={tracker}
      Actions={Actions}
    />
  );
}
function TrackBoxView(
  p: TrackProps & { box: RefObject<HTMLDivElement | null>; tracker: TrackBox; Actions: ComponentType<TrackProps> },
) {
  const { context: ctx, section, track, tracker, Actions } = p;
  // Focus inside it when it goes lands on the section's surface.
  useLanding(p.box);
  useLayoutEffect(() => {
    // The box is new: placed once now, then as the viewport or the layout moves.
    tracker.place();
    const [vp, layout] = [ctx.viewports.get(section), ctx.layout.get(section)];
    const subs = [
      vp ? later(vp.live.version, tracker.place) : undefined,
      layout ? later(layout.version, () => tracker.update()) : undefined,
    ];
    return () => off(subs);
  }, [ctx, section, tracker]);
  return (
    // Pointer and focus boundaries of a box that holds controls, not a control.
    // biome-ignore lint/a11y/noStaticElementInteractions: as above.
    // biome-ignore lint/a11y/noNoninteractiveElementInteractions: as above.
    <div
      ref={p.box}
      className="seqv-track-actions"
      onPointerEnter={() => tracker.hold(track)}
      onPointerLeave={() => tracker.hold(null)}
      onBlur={tracker.later}
    >
      <Actions context={ctx} section={section} track={track} row={tracker.row} />
    </div>
  );
}

interface RowsView {
  readonly top: number;
  readonly left: number;
  readonly width: number;
  readonly rows: readonly {
    readonly displayRow: number;
    readonly dataRow: number;
    readonly top: number;
    readonly height: number;
  }[];
}
/**
 * The `RowControl` slot beside each visible row of the section's rows layer, in the header's last pixels. The layout is
 * read in a subscription to the rows in view (its getters can build the layout and notify, which they must not do in a
 * render), edge-rate: when the rows in view change, not per scroll. The scroll moves the box by a custom property.
 */
function RowControls({ context: ctx, section, Control }: SectionProps & { Control: ComponentType<RowProps> }) {
  const [view, setView] = useState<RowsView | null>(null);
  const inner = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const rows = ctx.layout.visibleRows(section);
    const layout = ctx.layout.get(section);
    if (!rows || !layout) return;
    const sub = rows.version.subscribe(() => {
      const [r, l, w] = [rows.current, layout.current, ctx.layout.props(section).controlWidth];
      setView({
        top: l.scrolledTop,
        left: l.headerWidth - w,
        width: w,
        rows: Array.from({ length: r.count }, (_, i) => ({
          displayRow: r.displayRow[i],
          dataRow: r.dataRow[i],
          top: r.top[i],
          height: r.height[i],
        })),
      });
    });
    return () => sub.unsubscribe();
  }, [ctx, section]);
  const placed = view !== null;
  useLayoutEffect(() => {
    const vp = ctx.viewports.get(section);
    if (!placed || !vp) return;
    const at = { quarter: Number.NaN, page: Number.NaN };
    const sub = vp.live.version.subscribe(
      () => inner.current && placeScroll(inner.current, vp.live.current.scrollY, at),
    );
    return () => sub.unsubscribe();
  }, [ctx, section, placed]);
  if (!view) return null;
  return (
    <div className="seqv-row-controls" style={{ top: view.top, left: view.left, width: view.width }}>
      <div ref={inner}>
        {view.rows.map((r) => (
          <div
            key={`${r.displayRow}:${r.dataRow}`}
            className="seqv-row-control"
            style={{ top: r.top, height: r.height }}
          >
            <Control context={ctx} section={section} displayRow={r.displayRow} dataRow={r.dataRow} />
          </div>
        ))}
      </div>
    </div>
  );
}
