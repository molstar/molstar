/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from '../../dom/jsx-types.js';
import { type MaybeReactive, onCleanup } from '../../dom/scope.js';
import type { ShellModel } from '../../shell/shell.js';
import type { StatusModel } from '../../shell/status.js';
import { classes as canvasClasses } from './_generated/canvas.module.css.js';
import { classes } from './_generated/shell.module.css.js';
import { CanvasToolbar, OverlayHost, StatusLine } from './canvas-chrome.js';
import { DockView } from './dock-view.js';

export interface ShellViewProps {
  shell: ShellModel;
  status: StatusModel;
  /** The canvas (the Mol* viewport, or a placeholder in stories). */
  canvas: Child;
  /** Items floating over the canvas besides the toolbar and status line (snapshots, hover label, …). */
  overlays?: Child;
  /** The timeline at the bottom centre of the canvas (trajectory, animation). */
  timeline?: Child;
  /**
   * The theme of the chrome over the canvas, from the canvas background (light background: light glass), independent
   * of the UI theme: dark glass over a white canvas would be muddy.
   */
  canvasTheme?: MaybeReactive<'light' | 'dark'>;
  onOpenLog?: () => void;
}

/** The frame: docks around the canvas, chrome over it. */
export function ShellView(props: ShellViewProps) {
  const { shell, status } = props;
  const root = (
    <div class={classes.shell} data-mode={shell.mode}>
      <DockView shell={shell} dock={shell.dock('top')} />
      <DockView shell={shell} dock={shell.dock('left')} />
      <main class={classes.viewport}>
        <div class={classes.canvas}>{props.canvas}</div>
        <div class={classes.layer} data-msp-theme={props.canvasTheme ?? 'light'}>
          {props.overlays}
          <div class={canvasClasses.bottom}>
            <StatusLine status={status} onOpenLog={props.onOpenLog} />
            {props.timeline ? <div class={canvasClasses.timeline}>{props.timeline}</div> : null}
          </div>
          <OverlayHost shell={shell} />
          <CanvasToolbar shell={shell} />
        </div>
      </main>
      <DockView shell={shell} dock={shell.dock('right')} />
      <DockView shell={shell} dock={shell.dock('bottom')} />
    </div>
  ) as HTMLElement;
  // The shell decides compact or desktop by its own width, not the window's.
  const observer = new ResizeObserver((entries) => shell.setWidth(entries[0].contentRect.width));
  observer.observe(root);
  onCleanup(() => observer.disconnect());
  return root;
}
