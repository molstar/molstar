/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon, openContextMenu, type ReadonlyState } from '@molstar/prism-ui';
import { EntryKind, type SectionLayout } from '@molstar/seq-viewer';
import { getScope } from '@molstar/prism-ui/dom/scope';
import { switchMap } from 'rxjs';
import { type SequenceApiRunner, setTrackHidden, showAllTracks } from './commands.js';
import type { SequenceViewer, TrackColoring } from './model.js';
import { classes } from './_generated/track-gutter.module.css.js';

// The tracks' left gutter (sequence.md §8): a coloring icon per track row that has one (the action of the
// track's menu: pLDDT, secondary structure, the default coloring), a Prism DOM overlay over Seq*'s surface. The rows'
// tops and heights are read from Seq*'s section layout when the layout is rebuilt, never per frame or pointer move; the
// row the pointer is over is Seq*'s `interaction.hover` (its track, whatever part of the row), which fires only when
// what is under the pointer changes, and the DOM changes when the track differs. While tracks are hidden, a faint
// "N hidden" at the foot of the tracks opens a menu that shows them again (sequence.md §8); with none hidden it is not
// there.

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

/**
 * Mounts the gutter over `surface` (in `root`, which is positioned and lets the pointer through), for the rows of the
 * section the band shows (`section`); returns the unmount.
 */
export function mountTrackGutter(
  viewer: SequenceViewer,
  surface: HTMLElement,
  root: HTMLElement,
  section: ReadonlyState<string>,
  /** The API as the person: the hint's menu shows tracks through it (steps); none: on the model. */
  api?: SequenceApiRunner,
): () => void {
  const seq = viewer.seq;
  const gutter = (<div class={classes.gutter} />) as HTMLElement;
  root.append(gutter);

  // The hidden tracks' hint: its menu lists them (Show all first), opened under it.
  let closeMenu: (() => void) | undefined;
  const hint = (
    <button
      type="button"
      class={classes.hint}
      onClick={() => {
        const hidden = viewer.tracks.value.filter((t) => t.hidden);
        if (!hidden.length) return;
        const r = hint.getBoundingClientRect();
        closeMenu?.();
        closeMenu = openContextMenu(
          hint.closest<HTMLElement>('.msp-root') ?? document.body,
          { x: r.left, y: r.bottom },
          [
            { id: 'show-all', label: 'Show all tracks', icon: 'eye', run: () => void showAllTracks(viewer, api) },
            ...hidden.map((t) => ({
              id: `show:${t.ref}`,
              label: `Show ${t.label}`,
              run: () => void setTrackHidden(viewer, api, t.ref, false),
            })),
          ],
        );
      }}
    />
  ) as HTMLElement;
  /** The hint's text and place, written when they change. */
  let hintText = '';
  let hintTop = -1;
  const placeHint = (count: number, top: number) => {
    const text = count > 0 ? `${count} hidden` : '';
    if (text !== hintText) {
      hintText = text;
      hint.textContent = text;
      hint.setAttribute('aria-label', count === 1 ? '1 hidden track: show it' : `${count} hidden tracks: show them`);
      if (text) gutter.append(hint);
      else hint.remove();
    }
    if (text && top !== hintTop) {
      hintTop = top;
      hint.style.transform = `translateY(${top}px)`;
    }
  };

  const buttons = new Map<string, { el: HTMLElement; coloring: TrackColoring }>();
  let over: string | undefined;

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
    const layout = seq.layout.get(section.value)?.current;
    const view = a?.view;
    for (const { el } of buttons.values()) el.remove();
    buttons.clear();
    over = undefined;
    if (!layout || !view || viewer.style.value !== 'tracks') return placeHint(0, 0);
    const refOf = (t: number) => a.set.strings.get(view.resolved.tracks.ref[t]);
    // The hint sits under the last row shown.
    const all = gutterRows(layout, refOf, () => true);
    const last = all[all.length - 1];
    placeHint(viewer.tracks.value.filter((t) => t.hidden).length, last ? last.top + last.height : 0);
    const rows = gutterRows(layout, refOf, (ref) => viewer.trackColoring(ref) !== undefined);
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

  // Hidden UI does no UI updates (styling.md §5): while the band's scope is paused (a hidden dock), a change marks the
  // gutter stale and it is rebuilt once on resume, not on every layout version.
  let active = true;
  let stale = false;
  const update = () => {
    if (active) rebuild();
    else stale = true;
  };
  getScope()?.onActive((a) => {
    active = a;
    if (a && stale) {
      stale = false;
      rebuild();
    }
  });

  let stop: (() => void) | undefined;
  const follow = () => {
    stop?.();
    stop = undefined;
    const live = seq.layout.get(section.value);
    update();
    if (!live) return;
    const sub = live.version.subscribe(update);
    stop = () => sub.unsubscribe();
  };
  const sections = seq.viewports.sections.subscribe(follow);
  const style = viewer.style.subscribe(follow);
  const structure = viewer.structure.subscribe(follow);
  const shown = section.subscribe(follow);
  // Hiding or showing a track changes the hint (the layout rebuilds too, but the list may come after it).
  const tracks = viewer.tracks.subscribe(update);

  const hovers = section.pipe(switchMap((ref) => seq.interaction.on('hover', { section: ref })));
  const hovered = hovers.subscribe((e) =>
    hover(e.part.track !== null && buttons.has(e.part.track) ? e.part.track : undefined),
  );
  const onLeave = () => hover(undefined);
  // The icons are not inside the surface: the pointer going onto one is still over its row.
  const onSurfaceLeave = (e: PointerEvent) => {
    if (!gutter.contains(e.relatedTarget as Node | null)) onLeave();
  };
  const onGutterLeave = (e: PointerEvent) => {
    if (!surface.contains(e.relatedTarget as Node | null)) onLeave();
  };
  surface.addEventListener('pointerleave', onSurfaceLeave);
  gutter.addEventListener('pointerleave', onGutterLeave);

  return () => {
    stop?.();
    sections.unsubscribe();
    style.unsubscribe();
    structure.unsubscribe();
    shown.unsubscribe();
    tracks.unsubscribe();
    closeMenu?.();
    hovered.unsubscribe();
    surface.removeEventListener('pointerleave', onSurfaceLeave);
    gutter.removeEventListener('pointerleave', onGutterLeave);
    gutter.remove();
  };
}
