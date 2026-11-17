/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import {
  activity,
  bind,
  type Child,
  type MaybeReactive,
  observeResize,
  on,
  Show,
  setAttr,
  setStyle,
} from '../dom/index.js';
import { classes } from './_generated/stage.module.css.js';
import { EmptyStage } from './empty-stage.js';
import { HoverLabel } from './hover-label.js';
import type { StageModel } from './stage-model.js';

export interface StageProps {
  model: StageModel;
  /** The app's canvas (a WebGL canvas, an SVG, a document): it fills the stage; its role and name are the app's. */
  canvas: Child;
  /** A band over the canvas's top (the sequences): its backdrop the canvas's color, fading out below it. */
  top?: Child;
  /** The top left, below the band (a playback deck). */
  lead?: Child;
  /** The bottom left (the status line). */
  status?: Child;
  /** The bottom centre (what the person works on: chips). */
  chips?: Child;
  /** The bottom right (a `StageToolbar`). */
  corner?: Child;
  /**
   * Whether the corner's toolbar shows (stage §2; default always): the app decides, from what it knows (a pane's
   * status from its layout: shown while active or hovered, so several stages side by side show one toolbar, not one
   * each). Hidden, it fades out and takes no pointer.
   */
  toolbarShown?: MaybeReactive<boolean>;
  /**
   * What of its bottom edge something over it covers (px from each corner, stage §2): its bottom bar keeps clear (a
   * pane under a shell's shared chrome: its pane's insets, shell §2.1). The stage does not know what covers it.
   */
  insets?: MaybeReactive<{ readonly bottomLeft: number; readonly bottomRight: number }>;
  /** Anything else over the canvas (full size; only its content takes the pointer). */
  overlays?: Child;
  /** What shows on the empty stage (stage §4), the app's: default the Mol* glass splash; `() => null` for none. */
  emptyMark?: (ctx: { readonly compact: boolean }) => Child;
  /** What the pointer must be over for the hover label: the canvas and `data-msui-hover-source` (the band) by default. */
  hoverOver?: (target: Element) => boolean;
  class?: string;
}

const overCanvas = (t: Element) => !!t.closest?.('[data-msui-stage-canvas], [data-msui-hover-source]');

/** A pinch (ctrl + wheel, a trackpad's) zooms the page unless held: held where the stage zooms its own content. */
const holdPinch = (e: WheelEvent) => {
  if (e.ctrlKey) e.preventDefault();
};
const holdGesture = (e: Event) => e.preventDefault();

/**
 * The stage (stage §2): the app's canvas, and over it one layer of chrome in the canvas's ink (`StageModel.ink`: its
 * theme follows the canvas's color, not the UI's). The layer takes no pointer, only what is in its slots does: the band
 * on top (its faded part lets the pointer through), the empty mark, the lead below the band, the hover label, and the
 * bottom bar (status | chips | corner). It measures the band (`topInset`) and itself (`size`) into the model.
 */
export function Stage(props: StageProps) {
  const { model } = props;
  const canvasBox = (
    <div class={classes.canvas} data-msui-stage-canvas="">
      {props.canvas}
    </div>
  ) as HTMLElement;
  const bottom = (
    <div class={classes.bottom}>
      <div class={`${classes.slot} ${classes.status}`}>{props.status}</div>
      <div class={`${classes.slot} ${classes.chips}`}>{props.chips}</div>
      <div class={`${classes.slot} ${classes.corner}`}>
        <Show when={derive(model.hoverAt, (p) => p === 'corner')}>
          {() => <HoverLabel model={model.hover} area={() => root} over={props.hoverOver ?? overCanvas} />}
        </Show>
        {props.corner === undefined ? null : <div class={classes.tools}>{props.corner}</div>}
      </div>
    </div>
  ) as HTMLElement;
  const band =
    props.top === undefined
      ? undefined
      : (activity(model.band, () => (
          <div class={classes.top}>
            <div class={classes.band}>{props.top}</div>
          </div>
        )) as HTMLElement);
  const ink = (
    <div class={classes.ink}>
      {band}
      <EmptyStage model={model} mark={props.emptyMark} />
      {props.lead === undefined ? null : <div class={`${classes.slot} ${classes.lead}`}>{props.lead}</div>}
      <Show when={derive(model.hoverAt, (p) => p === 'cursor')}>
        {() => (
          <HoverLabel model={model.hover} area={() => root} over={props.hoverOver ?? overCanvas} placement="cursor" />
        )}
      </Show>
      {bottom}
      {props.overlays === undefined ? null : <div class={`${classes.slot} ${classes.overlays}`}>{props.overlays}</div>}
    </div>
  ) as HTMLElement;
  const root = (
    <div class={props.class ? `${classes.stage} ${props.class}` : classes.stage} data-msui-stage="">
      {canvasBox}
      {ink}
    </div>
  ) as HTMLElement;

  bind(model.ink, (v) => {
    setStyle(root, '--msui-stage-color', v.color);
    setAttr(ink, 'data-msui-theme', v.theme);
  });
  bind(model.quiet, (q) => setAttr(canvasBox, 'data-quiet', q));
  bind(model.compact, (c) => setAttr(root, 'data-compact', c));
  bind(props.insets ?? { bottomLeft: 0, bottomRight: 0 }, (i) => {
    setStyle(root, '--msui-stage-inset-left', `${i.bottomLeft}px`);
    setStyle(root, '--msui-stage-inset-right', `${i.bottomRight}px`);
  });
  // The toolbar shown or not is the app's (stage §2): it gives a value or a state, the stage follows it.
  bind(props.toolbarShown ?? true, (shown) => setAttr(root, 'data-toolbar', shown ? false : 'hidden'));
  bind(model.topInset, (px) => setStyle(root, '--msui-stage-top-inset', `${px}px`));
  if (band) bind(model.band, (shown) => setAttr(band, 'hidden', !shown));

  // Measured into the model: the stage's size, and the band's height (its content: the fade is drawn past it). One
  // observer each, its entries' sizes read (no layout forced).
  observeResize(root, (entries) => {
    const r = entries[entries.length - 1]?.contentRect;
    if (r) model.setSize(r.width, r.height);
  });
  if (band) {
    observeResize(band.firstElementChild ?? band, (entries) => {
      const r = entries[entries.length - 1]?.contentRect;
      if (r && model.band.value) model.setTopInset(r.height);
    });
  }

  // A pinch is held over the canvas and the bottom bar only (non-passive there, passive everywhere else: a wheel over
  // the band or a panel scrolls without waiting on a listener); Safari's gestures on the root.
  on<WheelEvent>(canvasBox, 'wheel', holdPinch, { passive: false });
  on<WheelEvent>(bottom, 'wheel', holdPinch, { passive: false });
  on(root, 'gesturestart', holdGesture);
  on(root, 'gesturechange', holdGesture);
  on(root, 'gestureend', holdGesture);
  return root;
}
