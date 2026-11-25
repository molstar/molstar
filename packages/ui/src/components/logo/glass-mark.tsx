/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { frames, observeIntersection, on, onActive, onCleanup } from '../../dom/index.js';
import { classes } from './_generated/glass-mark.module.css.js';
import { Glass, type GlassColors, type GlassTiming } from './glass.js';
import { onThemeChange } from './theme-change.js';
import { DefaultWordmarkLabels, WordmarkGeometry, type WordmarkLabels } from './wordmark.js';

/** The reduced-motion query, where there is one (followed: a change while the mark shows takes effect at once). */
function reducedMotionQuery(): MediaQueryList | undefined {
  return typeof matchMedia === 'function' ? matchMedia('(prefers-reduced-motion: reduce)') : undefined;
}

export interface GlassMarkProps {
  /** Whole pixels a unit. */
  scale: number;
  /** Room around the word for the pieces, in units (the canvas reaches this far past the word's box each way). */
  margin: number;
  /** The cut. */
  seed: number;
  /** How far the pieces travel, as a factor of the farthest (1). */
  travel: number;
  /** The assembly's timing (the glass's default without it). */
  timing?: GlassTiming;
  /** Whether it assembles as it shows (asked each time it shows); otherwise it shows whole. Default: always. */
  assemble?: () => boolean;
  class?: string;
  /** Its accessible name (`mark`: the word it draws). */
  labels?: Partial<WordmarkLabels>;
}

/**
 * The Mol* wordmark as glass (`glass.ts`), the one frame loop under `GlassWordmark` (a heading) and `GlassSplash` (a stage).
 * Each time it shows (its scope becomes active) it assembles from its pieces; they lean under the pointer, and a press
 * shatters it for a moment before it gathers again. It draws only while something moves, in the page's one frame tick
 * (`frames`): at rest, with a pointer resting over it, hidden (its scope paused) or off screen (`IntersectionObserver`)
 * it asks for no frames. Under reduced motion (followed as it changes), or with nothing to draw with, one still frame of
 * the word whole. Its colors (CSS) and its pixel ratio are read as it shows and again when the theme or the ratio
 * changes, never a frame.
 */
export function GlassMark(props: GlassMarkProps): Node {
  const scale = Math.max(1, Math.round(props.scale));
  const width = WordmarkGeometry.width * scale;
  const height = WordmarkGeometry.height * scale;
  const margin = props.margin * scale;
  const canvas = (
    <canvas
      class={classes.canvas}
      aria-hidden="true"
      style={{
        left: `${-margin}px`,
        top: `${-margin}px`,
        width: `${width + 2 * margin}px`,
        height: `${height + 2 * margin}px`,
      }}
    />
  ) as HTMLCanvasElement;
  // The colors are CSS (the logo's tokens, the subtle text: they follow the theme); read from these, not set and read
  // back on the canvas (a style write before every read).
  const probes = [1, 2, 3, 4, 0].map((n) => (<i class={classes.probe} data-n={String(n)} />) as HTMLElement);
  const root = (
    <span
      class={props.class ? `${classes.mark} ${props.class}` : classes.mark}
      style={{ width: `${width}px`, height: `${height}px` }}
      role="img"
      aria-label={props.labels?.mark ?? DefaultWordmarkLabels.mark}
    >
      {canvas}
      {probes}
    </span>
  ) as HTMLElement;

  const ctx = canvas.getContext('2d');
  const noDrawing = typeof Path2D !== 'function' || !ctx;
  const reduced = reducedMotionQuery();
  let still = noDrawing || !!reduced?.matches;
  let glass = new Glass(props.seed, props.travel);
  // The word's drawing origin in the canvas (its viewBox runs from -1).
  const origin = margin + scale;
  let dpr = 0;
  let colors: GlassColors = { logo: ['#000', '#000', '#000', '#000'], edge: '#000' };

  /** Sized for the screen (only when the pixel ratio changed) and colored from the theme where it shows. */
  const prepare = () => {
    if (!ctx) return;
    const ratio = Math.min(2, globalThis.devicePixelRatio || 1);
    if (ratio !== dpr) {
      dpr = ratio;
      canvas.width = Math.round((width + 2 * margin) * dpr);
      canvas.height = Math.round((height + 2 * margin) * dpr);
    }
    const read = probes.map((p) => getComputedStyle(p).color);
    colors = { logo: [read[0], read[1], read[2], read[3]], edge: read[4] };
  };

  const draw = () => {
    if (!ctx) return;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    glass.draw(ctx, origin, origin, scale, colors, dpr);
  };

  let last = 0;
  let active = true;
  /** On screen (true where there is no observer: jsdom, an old browser). */
  let onScreen = true;
  /** Frames wanted while it was off screen: asked for again when it comes back. */
  let wanted = false;
  /** Shown again: prepared and assembled at the next frame (once it is in the page). */
  let fresh = false;
  /** The theme or the pixel ratio changed: prepared again (colors, size) and drawn at the next frame. */
  let stale = false;
  const tick = (now: number) => {
    if (!active) return;
    if (!onScreen) {
      wanted = true;
      last = 0;
      return;
    }
    if (fresh) {
      fresh = false;
      prepare();
      if (!still && (props.assemble?.() ?? true)) glass.scatter(props.timing);
      last = 0;
    } else if (stale) {
      stale = false;
      prepare();
    }
    if (!still) {
      const dt = last ? Math.min(0.05, (now - last) / 1000) : 1 / 60;
      last = now;
      glass.step(dt);
    }
    draw();
    if (!still && !glass.resting) loop.request();
    else last = 0;
  };
  const loop = frames(tick);

  observeIntersection(root, (entries) => {
    const entry = entries[entries.length - 1];
    if (!entry) return;
    onScreen = entry.isIntersecting;
    if (!onScreen) {
      // Off screen: no frames; what was to come (its first frame, the pieces moving) comes when it is back.
      loop.cancel();
      if (fresh || !glass.resting) wanted = true;
    } else if (wanted) {
      wanted = false;
      loop.request();
    }
  });

  // Reduced motion turned on or off while it shows: on, the word whole and still at once (a new glass is whole); off,
  // it moves from the next pointer on.
  if (reduced && typeof reduced.addEventListener === 'function') {
    on(reduced, 'change', () => {
      const next = noDrawing || reduced.matches;
      if (next === still) return;
      still = next;
      if (still) glass = new Glass(props.seed, props.travel);
      last = 0;
      loop.request();
    });
  }

  // The theme or the pixel ratio changed (a theme switch, the window moved to another screen): read again as the next
  // frame draws (a resting glass is drawn once more, in its new colors).
  const restyle = () => {
    stale = true;
    loop.request();
  };
  if (!noDrawing) {
    onThemeChange(root, restyle);
    // The ratio's query is for the ratio now: made again at each change (no listener while it holds).
    let ratioQuery: MediaQueryList | undefined;
    const followRatio = () => {
      ratioQuery?.removeEventListener('change', onRatio);
      ratioQuery =
        typeof matchMedia === 'function'
          ? matchMedia(`(resolution: ${globalThis.devicePixelRatio || 1}dppx)`)
          : undefined;
      if (typeof ratioQuery?.addEventListener === 'function') ratioQuery.addEventListener('change', onRatio);
      else ratioQuery = undefined;
    };
    const onRatio = () => {
      followRatio();
      restyle();
    };
    followRatio();
    onCleanup(() => ratioQuery?.removeEventListener('change', onRatio));
  }

  if (!noDrawing) {
    // The pointer, in drawing units, from the box's corner as numbers: read when the pointer comes in and again after
    // anything scrolls or the window resizes (the box moved under a resting pointer), never per move; nothing is made
    // per move (no object, no rect). Its listeners do nothing while still (reduced motion may turn off).
    let boxLeft = 0;
    let boxTop = 0;
    let boxRead = false;
    const readBox = () => {
      const r = root.getBoundingClientRect();
      boxLeft = r.left;
      boxTop = r.top;
      boxRead = true;
    };
    const moved = () => {
      boxRead = false;
    };
    on(root, 'pointerenter', readBox);
    on(root.ownerDocument, 'scroll', moved, { capture: true, passive: true });
    if (root.ownerDocument.defaultView) on(root.ownerDocument.defaultView, 'resize', moved, { passive: true });
    on<PointerEvent>(root, 'pointermove', (e) => {
      if (still) return;
      if (!boxRead) readBox();
      glass.setPointer((e.clientX - boxLeft) / scale - 1, (e.clientY - boxTop) / scale - 1);
      loop.request();
    });
    on(root, 'pointerleave', () => {
      boxRead = false;
      if (still) return;
      glass.clearPointer();
      loop.request();
    });
    on<PointerEvent>(root, 'pointerdown', (e) => {
      if (still) return;
      // A press is the glass's: no text selection starts there (a double press selected the words around it). Not
      // stopped: what listens around it (a layer closing on a press outside) still hears it.
      e.preventDefault();
      if (!boxRead) readBox();
      glass.shatter((e.clientX - boxLeft) / scale - 1, (e.clientY - boxTop) / scale - 1);
      loop.request();
    });
  }

  // Each time it shows: prepared and assembled. Hidden: nothing (`frames` waits while the scope is paused).
  onActive((a) => {
    active = a;
    if (!a) {
      loop.cancel();
      last = 0;
      return;
    }
    fresh = true;
    loop.request();
  });
  onCleanup(() => {
    active = false;
  });
  return root;
}
