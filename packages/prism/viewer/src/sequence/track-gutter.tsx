/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon } from '@molstar/prism-ui';
import { EntryKind, type SectionLayout } from '@molstar/seq-viewer';
import { SectionRef } from './document.js';
import type { SequenceViewer, TrackColoring } from './model.js';
import { classes } from './_generated/track-gutter.module.css.js';

// The tracks' left gutter (prism-sequence.md §9, round 3): a coloring icon per track row that has one (the action of the
// track's menu: pLDDT, secondary structure, the default coloring), a Prism DOM overlay over Seq*'s surface. The rows'
// tops and heights are read from Seq*'s section layout when the layout is rebuilt, never per frame or pointer move; the
// pointer only compares the row it is over with the last, and the DOM changes when that differs.

/** A track row that has an icon: its place in the section (px from the surface's top). */
export interface GutterRow {
  readonly track: string;
  readonly top: number;
  readonly height: number;
}

/**
 * The rows of the section's tracks that `has` (by ref), from its layout. The pinned rows are at their surface y; the
 * scrolled ones start at `scrolledTop` (the tracks version does not scroll vertically: it is as tall as its content).
 */
export function gutterRows(
  layout: Pick<
    SectionLayout,
    'entryCount' | 'pinnedCount' | 'entryKind' | 'entryTrack' | 'entryTop' | 'entryHeight' | 'scrolledTop'
  >,
  refOf: (track: number) => string | undefined,
  has: (track: string) => boolean,
): GutterRow[] {
  const out: GutterRow[] = [];
  for (let i = 0; i < layout.entryCount; i++) {
    if (layout.entryKind[i] !== EntryKind.track) continue;
    const track = refOf(layout.entryTrack[i]);
    if (track === undefined || !has(track)) continue;
    const top = layout.entryTop[i] + (i < layout.pinnedCount ? 0 : layout.scrolledTop);
    out.push({ track, top, height: layout.entryHeight[i] });
  }
  return out;
}

/** Mounts the gutter over `surface` (in `root`, which is positioned and lets the pointer through); returns the unmount. */
export function mountTrackGutter(viewer: SequenceViewer, surface: HTMLElement, root: HTMLElement): () => void {
  const seq = viewer.seq;
  const gutter = (<div class={classes.gutter} />) as HTMLElement;
  root.append(gutter);

  let rows: GutterRow[] = [];
  const buttons = new Map<string, { el: HTMLElement; coloring: TrackColoring }>();
  let over: string | undefined;
  let rect: DOMRect | undefined;

  const refresh = () => {
    for (const { el, coloring } of buttons.values()) {
      const on = coloring.active();
      if (on !== el.hasAttribute('data-active')) el.toggleAttribute('data-active', on);
    }
  };
  const hover = (track: string | undefined) => {
    if (track === over) return;
    if (over !== undefined) buttons.get(over)?.el.removeAttribute('data-hover');
    over = track;
    if (track === undefined) return;
    buttons.get(track)?.el.setAttribute('data-hover', '');
    // The 3D view may have been colored elsewhere since: read it when the pointer comes onto another row.
    refresh();
  };

  /** The rows again, on the layout's rebuild (a resize, tracks coming and going, the style): edge rate. */
  const rebuild = () => {
    const a = seq.graph.assembled.value;
    const layout = seq.layout.get(SectionRef)?.current;
    const view = a?.view;
    gutter.replaceChildren();
    buttons.clear();
    rows = [];
    over = undefined;
    rect = undefined;
    if (!layout || !view || viewer.style.value !== 'tracks') return;
    rows = gutterRows(
      layout,
      (t) => a.set.strings.get(view.resolved.tracks.ref[t]),
      (ref) => viewer.trackColoring(ref) !== undefined,
    );
    for (const row of rows) {
      const coloring = viewer.trackColoring(row.track);
      const icon = coloring?.icon;
      if (!coloring || !icon) continue;
      const el = (
        <button
          type="button"
          class={classes.icon}
          aria-label={coloring.label}
          data-tip={coloring.label}
          data-tip-side="right"
          style={`transform: translateY(${row.top}px); height: ${row.height}px`}
          onClick={() => {
            coloring.run();
            // The action is asynchronous (one Mol* state update): the marker follows when it has run.
            setTimeout(refresh, 0);
            requestAnimationFrame(() => refresh());
          }}
        >
          <Icon name={icon} size={16} />
        </button>
      ) as HTMLElement;
      buttons.set(row.track, { el, coloring });
      gutter.append(el);
    }
    refresh();
  };

  let stop: (() => void) | undefined;
  const follow = () => {
    stop?.();
    stop = undefined;
    const live = seq.layout.get(SectionRef);
    rebuild();
    if (!live) return;
    const sub = live.version.subscribe(rebuild);
    stop = () => sub.unsubscribe();
  };
  const sections = seq.viewports.sections.subscribe(follow);
  const style = viewer.style.subscribe(follow);
  const structure = viewer.structure.subscribe(follow);

  const rowAt = (y: number) => {
    for (const r of rows) if (y >= r.top && y < r.top + r.height) return r.track;
    return undefined;
  };
  const onMove = (e: PointerEvent) => {
    rect ??= surface.getBoundingClientRect();
    hover(rowAt(e.clientY - rect.top));
  };
  const onLeave = () => {
    rect = undefined;
    hover(undefined);
  };
  // The icons are not inside the surface: the pointer going onto one is still over its row.
  const onSurfaceLeave = (e: PointerEvent) => {
    if (!gutter.contains(e.relatedTarget as Node | null)) onLeave();
  };
  const onGutterLeave = (e: PointerEvent) => {
    if (!surface.contains(e.relatedTarget as Node | null)) onLeave();
  };
  surface.addEventListener('pointermove', onMove);
  surface.addEventListener('pointerleave', onSurfaceLeave);
  gutter.addEventListener('pointerleave', onGutterLeave);

  return () => {
    stop?.();
    sections.unsubscribe();
    style.unsubscribe();
    structure.unsubscribe();
    surface.removeEventListener('pointermove', onMove);
    surface.removeEventListener('pointerleave', onSurfaceLeave);
    gutter.removeEventListener('pointerleave', onGutterLeave);
    gutter.remove();
  };
}
