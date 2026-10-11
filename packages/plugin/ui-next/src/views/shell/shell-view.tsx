/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from '../../dom/jsx-types.js';
import { onCleanup } from '../../dom/scope.js';
import type { ShellModel } from '../../shell/shell.js';
import type { StatusModel } from '../../shell/status.js';
import { classes } from './_generated/shell.module.css.js';
import { CanvasToolbar, OverlayHost, StatusLine } from './canvas-chrome.js';
import { DockView } from './dock-view.js';

export interface ShellViewProps {
  shell: ShellModel;
  status: StatusModel;
  /** The canvas (the Mol* viewport, or a placeholder in stories). */
  canvas: Child;
  /** Items floating over the canvas besides the toolbar and status line (playback, hover label, …). */
  overlays?: Child;
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
        <div class={classes.layer}>
          {props.overlays}
          <StatusLine status={status} onOpenLog={props.onOpenLog} />
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
