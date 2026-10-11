/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from Base UI's ScrollArea (https://github.com/mui/base-ui, MIT License, Copyright (c) 2019 Material-UI SAS):
 * native scrolling in a viewport with its own scrollbar hidden, an overlay thumb sized by the visible fraction and moved
 * by a transform, shown while hovered or scrolling, dragged with the pointer captured, a track press jumping there.
 */

import type { ElementProps } from '../../dom/jsx-types.js';
import { onCleanup } from '../../dom/scope.js';
import { classes } from './_generated/scroll-area.module.css.js';

/** The thumb never gets smaller than this (px). */
const MinThumb = 16;
/** Scrollbars stay shown this long after the last scroll (ms). */
const ScrollingMs = 500;

export type ScrollAreaProps = Omit<ElementProps<HTMLDivElement>, 'class'> & {
  /** The area: its size in the layout (flex, max-height), its surface and edges (background, radius, hairlines). */
  class?: string;
  /** The content inside: its layout (display, padding, gap). Role, aria, ref and events also go to the content. */
  contentClass?: string;
  /** Which way it scrolls: down (default), across, or both. */
  axis?: 'y' | 'x' | 'both';
  /** The scrolling element (to read or set its scroll position). */
  viewportRef?: (viewport: HTMLDivElement) => void;
};

/**
 * A scrolling area that looks the same everywhere: thin overlay scrollbars (no room taken, shown on hover and while
 * scrolling, hidden when nothing overflows) and edges that fade where there is more to scroll (reset.css). Every area
 * of the UI that scrolls uses it.
 */
export function ScrollArea(props: ScrollAreaProps) {
  const { class: rootClass, contentClass, axis = 'y', viewportRef, children, ...contentProps } = props;
  let root!: HTMLDivElement;
  let viewport!: HTMLDivElement;
  let content!: HTMLDivElement;
  const bars: { y?: Bar; x?: Bar } = {};

  // Geometry, recomputed on scroll and resize; sizes written only when they change.
  const update = () => {
    if (!viewport) return;
    if (bars.y) layoutBar(bars.y, viewport.clientHeight, viewport.scrollHeight, viewport.scrollTop);
    if (bars.x) layoutBar(bars.x, viewport.clientWidth, viewport.scrollWidth, viewport.scrollLeft);
  };
  const layoutBar = (bar: Bar, visible: number, total: number, position: number) => {
    const overflows = total - visible > 1;
    if (overflows !== bar.overflows) {
      bar.overflows = overflows;
      root.toggleAttribute(bar.axis === 'y' ? 'data-overflow-y' : 'data-overflow-x', overflows);
    }
    if (!overflows) return;
    const track = bar.axis === 'y' ? bar.track.clientHeight : bar.track.clientWidth;
    const size = Math.max(MinThumb, Math.round((track * visible) / total));
    if (size !== bar.size) {
      bar.size = size;
      if (bar.axis === 'y') bar.thumb.style.height = `${size}px`;
      else bar.thumb.style.width = `${size}px`;
    }
    bar.travel = track - size;
    bar.range = total - visible;
    const offset = Math.round((bar.travel * Math.min(Math.max(position, 0), bar.range)) / bar.range);
    if (offset !== bar.offset) {
      bar.offset = offset;
      bar.thumb.style.transform = bar.axis === 'y' ? `translateY(${offset}px)` : `translateX(${offset}px)`;
    }
  };

  // Shown while scrolling, for a moment after.
  let scrollingTimer: ReturnType<typeof setTimeout> | undefined;
  const stopScrolling = () => {
    scrollingTimer = undefined;
    root.removeAttribute('data-scrolling');
  };
  const onScroll = () => {
    update();
    if (scrollingTimer === undefined) root.setAttribute('data-scrolling', '');
    else clearTimeout(scrollingTimer);
    scrollingTimer = setTimeout(stopScrolling, ScrollingMs);
  };

  const scrollTo = (bar: Bar, value: number) => {
    if (bar.axis === 'y') viewport.scrollTop = value;
    else viewport.scrollLeft = value;
  };
  const scrollOf = (bar: Bar) => (bar.axis === 'y' ? viewport.scrollTop : viewport.scrollLeft);
  const along = (bar: Bar, e: PointerEvent) => (bar.axis === 'y' ? e.clientY : e.clientX);

  // A press on the thumb drags it; on the track it jumps there (the thumb centred on the pointer) and drags on.
  const onPointerDown = (bar: Bar) => (e: PointerEvent) => {
    if (e.button !== 0 || !bar.overflows) return;
    e.preventDefault();
    if (e.target !== bar.thumb && bar.travel > 0) {
      const r = bar.track.getBoundingClientRect();
      const at = along(bar, e) - (bar.axis === 'y' ? r.top : r.left) - bar.size / 2;
      scrollTo(bar, (Math.min(Math.max(at, 0), bar.travel) / bar.travel) * bar.range);
    }
    const start = along(bar, e);
    const startScroll = scrollOf(bar);
    bar.track.setPointerCapture(e.pointerId);
    bar.track.setAttribute('data-dragging', '');
    const move = (ev: PointerEvent) => {
      // A move without the button held: the release went missing (capture dropped); end here.
      if ((ev.buttons & 1) === 0) return end();
      if (bar.travel > 0) scrollTo(bar, startScroll + ((along(bar, ev) - start) / bar.travel) * bar.range);
    };
    const end = () => {
      bar.track.removeEventListener('pointermove', move);
      bar.track.removeEventListener('pointerup', end);
      bar.track.removeEventListener('pointercancel', end);
      bar.track.removeAttribute('data-dragging');
    };
    bar.track.addEventListener('pointermove', move);
    bar.track.addEventListener('pointerup', end);
    bar.track.addEventListener('pointercancel', end);
  };

  const scrollbar = (which: 'y' | 'x') => {
    const thumb = (<div class={classes.thumb} />) as HTMLDivElement;
    const track = (
      <div
        class={classes.scrollbar}
        data-orientation={which === 'y' ? 'vertical' : 'horizontal'}
        aria-hidden="true"
        // Like a native scrollbar, a press does not take focus.
        onMouseDown={(e) => e.preventDefault()}
      >
        {thumb}
      </div>
    ) as HTMLDivElement;
    const bar: Bar = { axis: which, track, thumb, overflows: false, size: -1, offset: -1, travel: 0, range: 0 };
    track.addEventListener('pointerdown', onPointerDown(bar));
    bars[which] = bar;
    return track;
  };

  root = (
    <div class={`${classes.root} ${rootClass ?? ''}`} data-axis={axis}>
      <div
        class={classes.viewport}
        ref={(el) => {
          viewport = el;
          viewportRef?.(el);
        }}
        onScroll={onScroll}
      >
        <div
          {...contentProps}
          class={`${classes.content} ${contentClass ?? ''}`}
          ref={(el: HTMLDivElement) => {
            content = el;
            contentProps.ref?.(el);
          }}
        >
          {children}
        </div>
      </div>
      {axis !== 'x' ? scrollbar('y') : null}
      {axis !== 'y' ? scrollbar('x') : null}
    </div>
  ) as HTMLDivElement;

  // The viewport (the area resized) and the content (items added, opened, removed) change what overflows.
  if (typeof ResizeObserver !== 'undefined') {
    const observer = new ResizeObserver(update);
    observer.observe(viewport);
    observer.observe(content);
    onCleanup(() => observer.disconnect());
  }
  onCleanup(() => {
    if (scrollingTimer !== undefined) clearTimeout(scrollingTimer);
  });
  queueMicrotask(update);
  return root;
}

interface Bar {
  axis: 'y' | 'x';
  track: HTMLDivElement;
  thumb: HTMLDivElement;
  overflows: boolean;
  /** Thumb size, offset, how far it can move (px), and how far the content can scroll (px). */
  size: number;
  offset: number;
  travel: number;
  range: number;
}
