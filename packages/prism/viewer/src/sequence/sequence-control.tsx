/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import { BandFadePx, bind, type DockId, onCleanup, openContextMenu, type ShellModel } from '@molstar/prism-ui';
import type { HoverLabel } from '@molstar/seq-viewer';
import { mountSection } from '@molstar/seq-viewer/canvas';
import { combineLatest, of, switchMap } from 'rxjs';
import { SectionRef } from './document.js';
import { mountTrackGutter } from './track-gutter.js';
import { linkSequence, type SequenceView3D } from './link.js';
import type { SequenceViewer } from './model.js';
import { classes } from './_generated/sequence-control.module.css.js';

// The sequence control (prism-sequence.md §4): Seq*'s section surface in the canvas band, the link with the 3D view
// while it is shown, the Seq* hover label for the shell's one cursor-following label, and the fades where the text
// scrolls. It has no controls of its own: the style, the colors and the legend are the canvas stack's Sequence menu
// (sequence/menu.ts, commands). Seq* draws the letters, hover, selection and focus and handles the pointer.

export interface SequenceControlProps {
  viewer: SequenceViewer;
  /**
   * The 3D view of a structure, for the link (hover and selection both ways, focus); none: the view stands alone.
   * Over the real plugin it is `pluginView3D(plugin, structure)` (sources/).
   */
  view3d?: (structure: Structure) => SequenceView3D | undefined;
  /** The shell and its dock, for the height rules (§4): tracks follow their content, the text has a fixed height. */
  shell?: ShellModel;
  dock?: DockId;
  /** The text version's lines' height, px (the band is this and its air and fade). */
  textHeight?: number;
  /** The Seq* hover label as lines for the shell's hover label (`[]`: none), once a frame and only when it changes. */
  onHover?: (lines: readonly string[]) => void;
}

/** The keyboard cursor's place is said this long after the last move, ms (not every step of a held arrow). */
const CursorSayMs = 300;

/** The text version's lines' height: two or three lines of letters; the rest scrolls. */
export const TextHeight = 80;
/** The air above the lines, px (`--msp-space-3`). */
const TopPadPx = 8;
/** What the band takes besides its lines: the air above them and the faded bottom (`--msp-band-fade`). */
const BandExtra = TopPadPx + BandFadePx;
/** The scrollbar stays this long after the last scroll step, ms. */
const ScrollingMs = 500;
/** Px the text's edges fade over where there is more to scroll. */
const FadePx = 18;

const NoLines: readonly string[] = [];

/** A Seq* hover label as the shell's lines: the title, then `key: value` per line. */
export function labelLines(label: HoverLabel | null): readonly string[] {
  if (!label) return NoLines;
  const lines = [label.title];
  for (const l of label.lines) lines.push(l.value ? `${l.key}: ${l.value}` : l.key);
  return lines;
}

const sameLines = (a: readonly string[], b: readonly string[]) =>
  a.length === b.length && a.every((x, i) => x === b[i]);

export function SequenceControl(props: SequenceControlProps) {
  const { viewer } = props;
  const seq = viewer.seq;

  const surface = (<div class={classes.surface} />) as HTMLElement;
  const section = mountSection(seq, SectionRef, surface);
  onCleanup(() => section.dispose());

  // The fonts: the host's tokens, read once at the start of a frame (the control is in the document by then).
  const fontFrame = requestAnimationFrame(() => {
    const css = getComputedStyle(surface);
    const mono = css.getPropertyValue('--msp-font-mono').trim();
    const ui = css.getPropertyValue('--msp-font-sans').trim();
    const size = Number.parseFloat(css.getPropertyValue('--msp-font-size-sm'));
    if (mono && ui) viewer.setFonts({ mono, ui, ...(size > 0 ? { size } : {}) });
  });
  onCleanup(() => cancelAnimationFrame(fontFrame));

  // The link while this control is rendered (link.ts): subscribed on mount, dropped on cleanup.
  if (props.view3d) onCleanup(linkSequence(viewer, props.view3d));

  // The track menus (prism-sequence.md §4): the browser's context-menu request on a track (a right-click, Ctrl+click on
  // macOS, the menu key or Shift+F10 at Seq*'s keyboard cursor) is Seq*'s `interaction.menu`, with the track's ref and
  // the place; the model says what the track offers, Prism's menu opens there and runs it on the 3D view. The browser's
  // own menu is suppressed only over a track that offers one.
  {
    let close: (() => void) | undefined;
    const menus = seq.interaction.on('menu', { section: SectionRef, part: ['data', 'header', 'gutter'] });
    const asked = menus.subscribe((e) => {
      const actions = e.part.track === null ? [] : viewer.trackActions(e.part.track);
      if (actions.length === 0) return;
      e.native.preventDefault();
      close?.();
      const root = surface.closest<HTMLElement>('.msp-root') ?? document.body;
      close = openContextMenu(
        root,
        { x: e.place.clientX, y: e.place.clientY },
        actions.map((a) => ({ id: a.id, label: a.label, icon: a.icon, run: a.run })),
      );
    });
    onCleanup(() => {
      asked.unsubscribe();
      close?.();
    });
  }

  // The hover label: the latest label of a frame, handed on when its lines differ from the last handed on.
  if (props.onHover) {
    const handOn = props.onHover;
    let pending: HoverLabel | null = null;
    let last = NoLines;
    let frame = 0;
    const flush = () => {
      frame = 0;
      const lines = labelLines(pending);
      if (sameLines(lines, last)) return;
      last = lines;
      handOn(lines);
    };
    const sub = seq.labels.current.subscribe((label) => {
      pending = label;
      if (!frame) frame = requestAnimationFrame(flush);
    });
    // Leaving the surface ends the label at once (Seq* clears its hover on the same event).
    const leave = () => {
      pending = null;
      if (!frame) frame = requestAnimationFrame(flush);
    };
    surface.addEventListener('pointerleave', leave);
    onCleanup(() => {
      sub.unsubscribe();
      surface.removeEventListener('pointerleave', leave);
      if (frame) cancelAnimationFrame(frame);
      if (last.length > 0) handOn(NoLines);
    });
  }

  // The fades: the text's rows scroll inside Seq*'s canvas, so the surface's mask is driven from its viewport state
  // (CSS cannot see that scroll). A side fades where there is more to scroll that way; the custom properties are written
  // when a side turns on or off, never per scroll step.
  {
    let top = false;
    let bottom = false;
    let stopVp: (() => void) | undefined;
    const fade = (t: boolean, b: boolean) => {
      if (t === top && b === bottom) return;
      if (t !== top) surface.style.setProperty('--seq-fade-top', t ? `${FadePx}px` : '0px');
      if (b !== bottom) surface.style.setProperty('--seq-fade-bottom', b ? `${FadePx}px` : '0px');
      [top, bottom] = [t, b];
      if (t || b) surface.dataset.fade = '';
      else delete surface.dataset.fade;
    };
    const follow = () => {
      stopVp?.();
      stopVp = undefined;
      const vp = seq.viewports.get(SectionRef);
      const layout = seq.layout.get(SectionRef);
      if (!vp || !layout) return fade(false, false);
      const sub = vp.live.version.subscribe(() => {
        const s = vp.live.current;
        // Only the text scrolls vertically here; the visible rows are the viewport less the pinned part above them.
        const text = viewer.style.value === 'compact';
        const rest = s.contentHeight - s.scrollY - (s.height - layout.current.scrolledTop);
        fade(text && s.scrollY > 0.5, text && rest > 0.5);
      });
      stopVp = () => sub.unsubscribe();
    };
    const sections = seq.viewports.sections.subscribe(follow);
    const style = viewer.style.subscribe(follow);
    onCleanup(() => {
      sections.unsubscribe();
      style.unsubscribe();
      stopVp?.();
    });
  }

  // The scrollbar: shown while the pointer is over the band, the text scrolls (for a moment after) or a drag is on; the
  // attribute is written when that changes, never per move or per scroll step.
  {
    let hover = false;
    let dragging = false;
    let scrolling = false;
    let lastScroll = 0;
    let timer: ReturnType<typeof setTimeout> | undefined;
    const show = () => {
      const on = hover || dragging || scrolling;
      if (on !== surface.hasAttribute('data-bars')) surface.toggleAttribute('data-bars', on);
    };
    const settle = () => {
      timer = undefined;
      const rest = ScrollingMs - (performance.now() - lastScroll);
      if (rest > 0) timer = setTimeout(settle, rest);
      else {
        scrolling = false;
        show();
      }
    };
    const enter = () => ((hover = true), show());
    const leave = () => ((hover = false), show());
    surface.addEventListener('pointerenter', enter);
    surface.addEventListener('pointerleave', leave);
    const drag = seq.ui.drag.subscribe((d) => ((dragging = d !== null), show()));
    // The scroll: the viewport's `scrollY` moving (the section may be replaced, so it is followed by the sections).
    let lastY = 0;
    let stopScroll: (() => void) | undefined;
    const follow = () => {
      stopScroll?.();
      stopScroll = undefined;
      const vp = seq.viewports.get(SectionRef);
      if (!vp) return;
      const sub = vp.live.version.subscribe(() => {
        const y = vp.live.current.scrollY;
        if (y === lastY) return;
        lastY = y;
        lastScroll = performance.now();
        if (!scrolling) {
          scrolling = true;
          show();
        }
        if (timer === undefined) timer = setTimeout(settle, ScrollingMs);
      });
      stopScroll = () => sub.unsubscribe();
    };
    const sections = seq.viewports.sections.subscribe(follow);
    onCleanup(() => {
      surface.removeEventListener('pointerenter', enter);
      surface.removeEventListener('pointerleave', leave);
      drag.unsubscribe();
      sections.unsubscribe();
      stopScroll?.();
      if (timer !== undefined) clearTimeout(timer);
    });
  }

  // The height rules: tracks follow Seq*'s content height, the text has its own; the dock is restored when the control
  // goes. Written only when the size changes.
  if (props.shell) {
    const { shell } = props;
    const dockId = props.dock ?? 'top';
    const before = shell.dock(dockId).state.value.size;
    const text = props.textHeight ?? TextHeight;
    onCleanup(() => shell.setSize(dockId, before));
    const content = combineLatest([viewer.style, viewer.structure]).pipe(
      switchMap(([style]) => (style === 'tracks' ? (seq.layout.contentHeight(SectionRef) ?? of(0)) : of(0))),
    );
    bind(combineLatest([viewer.style, content]), ([style, height]) => {
      shell.setSize(dockId, (style === 'tracks' && height > 0 ? height : text) + BandExtra);
    });
  }

  // The keyboard cursor, said (README §3.7): where Seq*'s cursor is ("A 23 LYS, …", its hover label while the cursor
  // holds the hover), once the arrows rest, in a polite region of the control's own; nothing while the pointer hovers.
  const cursorSaid = (<span class={classes.said} aria-live="polite" />) as HTMLElement;
  {
    let timer: ReturnType<typeof setTimeout> | undefined;
    const say = () => {
      timer = undefined;
      const text =
        seq.ui.keyCursor.value && seq.hover.live.current.owner === 'keyboard'
          ? labelLines(seq.labels.current.value).join(', ')
          : '';
      if (cursorSaid.textContent !== text) cursorSaid.textContent = text;
    };
    const later = () => {
      if (timer !== undefined) clearTimeout(timer);
      timer = setTimeout(say, CursorSayMs);
    };
    const moved = seq.ui.keyCursor.subscribe(later);
    // The labels follow the pointer too: only the cursor's moves are timed (nothing per pointer move).
    const labelled = seq.labels.current.subscribe(() => {
      if (seq.ui.keyCursor.value || cursorSaid.textContent) later();
    });
    onCleanup(() => {
      moved.unsubscribe();
      labelled.unsubscribe();
      if (timer !== undefined) clearTimeout(timer);
    });
  }

  const root = (
    <div class={classes.root}>
      {surface}
      {cursorSaid}
    </div>
  ) as HTMLElement;
  // The coloring icons of the tracks (round 3), where the host offers colorings.
  onCleanup(mountTrackGutter(viewer, surface, root));
  return root;
}
