/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { type MaybeReactive, onCleanup } from '../../dom/scope.js';
import type { ShellModel } from '../../shell/shell.js';
import type { StatusModel } from '../../shell/status.js';
import { classes as canvasClasses } from './_generated/canvas.module.css.js';
import { classes } from './_generated/shell.module.css.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import { CanvasStack, CanvasToolbar, ContextChips, OverlayHost, StatusLine } from './canvas-chrome.js';
import { DockView, FloatingDock } from './dock-view.js';

/** Widths (px) of what floats on the right: the rail (tab strip, toolbar), the gap to the card, the edge margin. */
const RailWidth = 44;
const Gap = 8;
const Edge = 12;

export interface ShellViewProps {
  shell: ShellModel;
  status: StatusModel;
  /** The canvas (the Mol* viewport, or a placeholder in stories). */
  canvas: Child;
  /** Items floating over the canvas besides the toolbar and status line (snapshots, hover label, …). */
  overlays?: Child;
  /** Playback at the top left of the canvas (the snapshot and trajectory card). */
  playback?: Child;
  /**
   * The theme of the chrome over the canvas, from the canvas background (light background: light glass), independent
   * of the UI theme: dark glass over a white canvas would be muddy.
   */
  canvasTheme?: MaybeReactive<'light' | 'dark'>;
  /**
   * The canvas background (a CSS color). The chrome over the canvas is tinted from it instead of the UI's own
   * background, so cards read as frosted canvas rather than slabs of panel color.
   */
  canvasColor?: MaybeReactive<string>;
  /** The top of the canvas stack, above the canvas toolbar items: the interact / select mode buttons. */
  canvasLead?: Child;
  /** Beside the canvas stack, outside its frame: what opens from it (the picking level, next to Select). */
  canvasAside?: Child;
  onOpenLog?: () => void;
}

/** The frame: docks around the canvas, chrome over it. */
export function ShellView(props: ShellViewProps) {
  const { shell, status } = props;
  const right = shell.dock('right');
  const top = shell.dock('top');
  // The right dock floats over the canvas unless it is docked, or the shell is compact (then it is a bottom sheet).
  const floating = combineLatest([right.state, shell.mode]).pipe(
    map(([s, mode]) => s.floating && mode !== 'compact'),
    distinctUntilChanged(),
  );
  const floatingTop = combineLatest([top.state, shell.mode]).pipe(
    map(([s, mode]) => s.floating && mode !== 'compact'),
    distinctUntilChanged(),
  );
  // What floats at the top (the sequence card) takes this much: the snapshots and the canvas stack go below.
  const topInset = combineLatest([top.state, top.active, floatingTop]).pipe(
    map(([s, active, f]) => (f && active ? Edge + s.size + Gap : Edge)),
    distinctUntilChanged(),
    map((px) => `${px}px`),
  );
  // What floats on the right takes this much of the canvas: the bottom bar keeps clear of it.
  const rightInset = combineLatest([right.state, right.active, floating]).pipe(
    map(([s, active, f]) => RailWidth + Gap + (f && active ? s.size + Gap : 0)),
    distinctUntilChanged(),
    map((px) => `${px}px`),
  );
  // An open floating right panel: overlays then open in its place instead of stacking glass on glass.
  const rightPanelOpen = combineLatest([right.active, floating]).pipe(
    map(([active, f]) => (f && active ? 'open' : undefined)),
    distinctUntilChanged(),
  );
  const rightPanelWidth = right.state.pipe(map((s) => `${s.size}px`));
  const root = (
    <div class={classes.shell} data-mode={shell.mode}>
      <Show when={floatingTop.pipe(map((f) => !f))}>{() => <DockView shell={shell} dock={top} />}</Show>
      <DockView shell={shell} dock={shell.dock('left')} />
      <main class={classes.viewport}>
        <div class={classes.canvas}>{props.canvas}</div>
        <div
          class={classes.layer}
          data-msp-theme={props.canvasTheme ?? 'light'}
          data-right-panel={rightPanelOpen}
          data-overlay={shell.overlay.pipe(map((o) => (o ? 'open' : undefined)))}
          style={{
            // Every derived shade (glass, controls, hover, muted text) follows from the base color.
            '--msp-color-bg': props.canvasColor,
            '--msp-canvas-right-inset': rightInset,
            '--msp-canvas-top-inset': topInset,
            // Clear of the rail only (what sits beside it whatever the panel does: the compact stack).
            '--msp-canvas-rail-inset': `${RailWidth + Gap}px`,
            '--msp-right-panel-width': rightPanelWidth,
            // Where cards beside the rail end (the floating panel's right edge).
            '--msp-canvas-panel-right': `${RailWidth + Gap}px`,
          }}
        >
          <Show when={floatingTop}>{() => <FloatingDock shell={shell} dock={top} />}</Show>
          <div class={canvasClasses.lead}>{props.playback}</div>
          {props.overlays}
          <div class={canvasClasses.bottom}>
            <div class={canvasClasses.bottomGrid}>
              <StatusLine status={status} onOpenLog={props.onOpenLog} />
              <ContextChips status={status} />
            </div>
          </div>
          <CanvasStack shell={shell} lead={props.canvasLead} aside={props.canvasAside} />
          <OverlayHost shell={shell} />
          <Show
            when={floating}
            fallback={() => (
              <div class={`${canvasClasses.toolbarSlot} ${canvasClasses.frame}`}>
                <CanvasToolbar shell={shell} />
              </div>
            )}
          >
            {() => (
              <FloatingDock
                shell={shell}
                dock={right}
                railItems={<CanvasToolbar shell={shell} part="rest" />}
                railFoot={<CanvasToolbar shell={shell} part="lead" />}
              />
            )}
          </Show>
        </div>
      </main>
      <Show when={floating.pipe(map((f) => !f))}>{() => <DockView shell={shell} dock={right} />}</Show>
      <DockView shell={shell} dock={shell.dock('bottom')} />
    </div>
  ) as HTMLElement;
  // Escape closes the open right panel, after anything nearer has had it: handlers on the document (an overlay card, an
  // open list) run first and mark it handled. Not while typing in a field, where Escape belongs to the field.
  const onKeyDown = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented || shell.overlay.value || isEditable(e.target)) return;
    if (!right.state.value.collapsed && right.active.value !== undefined) {
      shell.toggle('right');
      e.preventDefault();
    }
  };
  window.addEventListener('keydown', onKeyDown);
  onCleanup(() => window.removeEventListener('keydown', onKeyDown));
  // Tooltips for the icon buttons of the rail and the canvas stack.
  onCleanup(installTooltips(root));
  // The shell decides compact or desktop by its own width, not the window's.
  const observer = new ResizeObserver((entries) => shell.setWidth(entries[0].contentRect.width));
  observer.observe(root);
  onCleanup(() => observer.disconnect());
  return root;
}

function isEditable(target: EventTarget | null): boolean {
  if (!(target instanceof HTMLElement)) return false;
  if (target.isContentEditable) return true;
  if (target instanceof HTMLTextAreaElement || target instanceof HTMLSelectElement) return true;
  return target instanceof HTMLInputElement && !['checkbox', 'radio', 'range', 'button', 'color'].includes(target.type);
}
