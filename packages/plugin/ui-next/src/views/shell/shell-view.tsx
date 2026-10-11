/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable } from 'rxjs';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { type MaybeReactive, onCleanup } from '../../dom/scope.js';
import type { ShellModel } from '../../shell/shell.js';
import type { StatusModel } from '../../shell/status.js';
import { classes as canvasClasses } from './_generated/canvas.module.css.js';
import { classes } from './_generated/shell.module.css.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import {
  BusyOverlay,
  CanvasStack,
  CanvasToolbar,
  ContextChips,
  OverlayHost,
  RailMenu,
  StatusLine,
} from './canvas-chrome.js';
import { DockView, FloatingDock } from './dock-view.js';
import { EmptyCanvas } from './empty-canvas.js';
import { FileDropOverlay } from './file-drop.js';

/**
 * The rail's width, and the one margin of the canvas frame (px, `--msp-canvas-edge`): every floating thing keeps it
 * from the canvas edges, the docked rail, the open panel, and the sequence card.
 */
const RailWidth = 44;
const Edge = 8;

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
   * The work in hand at the bottom right of the canvas, before the status chips (the selection bar). A render
   * function: the compact layout docks it elsewhere, and each place renders its own.
   */
  context?: () => Child;
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
  /** The command menu while it is open: over the whole shell. */
  commandMenu?: Child;
  /** The top of the canvas stack, above the canvas toolbar items: the interact / select mode buttons. */
  canvasLead?: Child;
  /** Beside the canvas stack, outside its frame: what opens from it (the picking level, next to Select). */
  canvasAside?: Child;
  onOpenLog?: () => void;
  /** Nothing loaded: a while after, the Mol* wordmark shows in the middle of the canvas. */
  empty?: Observable<boolean>;
  /**
   * Files dropped anywhere on the plugin (to the plugin's drag-and-drop manager). Given, dragging files in shows a drop
   * overlay over everything; the host leaves it out when the spec turns the overlay off.
   */
  onDropFiles?: (files: readonly File[]) => void;
  /** File formats the drop reads, for the overlay's hint. */
  dropFormats?: readonly string[];
}

/** The frame: docks around the canvas, chrome over it. */
export function ShellView(props: ShellViewProps) {
  const { shell, status } = props;
  const right = shell.dock('right');
  const compact = shell.mode.pipe(
    map((m) => m === 'compact'),
    distinctUntilChanged(),
  );
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
    map(([s, active, f]) => (f && active ? Edge + s.size + Edge : Edge)),
    distinctUntilChanged(),
    map((px) => `${px}px`),
  );
  // Pinned: the canvas ends where the rail, or the open panel, begins (they sit beside it, not over it).
  const pinned = combineLatest([right.state, floating]).pipe(
    map(([s, f]) => f && s.pinned),
    distinctUntilChanged(),
  );
  // An overlay card open in the panel's place counts as the panel.
  const pinnedWidth = combineLatest([right.state, right.active, shell.overlay, pinned]).pipe(
    map(([s, active, overlay, p]) => (p ? RailWidth + (active || overlay ? s.size : 0) : 0)),
    distinctUntilChanged(),
  );
  // What floats on the right takes this much of the canvas: the bottom bar keeps clear of it.
  const rightInset = combineLatest([right.state, right.active, floating, pinned]).pipe(
    // The docked rail is the canvas's right edge, and the open panel is too: what sits beside either keeps the frame
    // margin. Pinned, both are outside the canvas: only the margin is left.
    map(([s, active, f, p]) => (p ? Edge : RailWidth + Edge + (f && active ? s.size + Edge : 0))),
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
      <main
        class={classes.viewport}
        style={{
          '--msp-canvas-edge': `${Edge}px`,
          '--msp-canvas-right-inset': rightInset,
          '--msp-canvas-top-inset': topInset,
          // Clear of the rail only (what sits beside it whatever the panel does: the compact stack).
          '--msp-canvas-rail-inset': `${RailWidth + Edge}px`,
          '--msp-right-panel-width': rightPanelWidth,
          // Where cards beside the rail end (the floating panel's right edge).
          '--msp-canvas-panel-right': `${RailWidth + Edge}px`,
          // How much of the right the pinned rail and panel take from the canvas (0 when floating over it).
          '--msp-canvas-pinned': pinnedWidth.pipe(map((px) => `${px}px`)),
        }}
      >
        <div class={classes.canvas}>{props.canvas}</div>
        {/* Chrome that is part of the view: themed by the canvas background, tinted from its color. */}
        <div
          class={`${classes.layer} ${classes.canvasLayer}`}
          data-msp-theme={props.canvasTheme ?? 'light'}
          style={{
            // Every derived shade (glass, controls, hover, muted text) follows from the base color.
            '--msp-color-bg': props.canvasColor,
          }}
        >
          <Show when={floatingTop}>{() => <FloatingDock shell={shell} dock={top} />}</Show>
          {props.empty ? <EmptyCanvas empty={props.empty} /> : null}
          <div class={canvasClasses.lead}>{props.playback}</div>
          {props.overlays}
          <Show when={compact.pipe(map((c) => !c))}>
            {() => (
              <div class={canvasClasses.bottom}>
                <div class={canvasClasses.bottomGrid}>
                  <StatusLine status={status} onOpenLog={props.onOpenLog} />
                  <ContextChips status={status} lead={props.context?.()} />
                </div>
              </div>
            )}
          </Show>
          <CanvasStack shell={shell} lead={props.canvasLead} aside={props.canvasAside} />
        </div>
        {/* The application's own chrome over the canvas (the rail, the panel, the cards): the UI theme. */}
        <div
          class={classes.layer}
          data-right-panel={rightPanelOpen}
          data-pinned={pinned.pipe(map((p) => (p ? '' : undefined)))}
          data-overlay={shell.overlay.pipe(map((o) => (o ? 'open' : undefined)))}
        >
          <OverlayHost shell={shell} />
          <Show
            when={floating}
            fallback={() => (
              // Docked beside the canvas (not compact: there the rail's actions are under More).
              <Show when={compact.pipe(map((c) => !c))}>
                {() => (
                  <div class={`${canvasClasses.toolbarSlot} ${canvasClasses.frame}`}>
                    <CanvasToolbar shell={shell} />
                  </div>
                )}
              </Show>
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
      {/* Compact: progress and the selection docked under the canvas, not over it. */}
      <Show when={compact}>
        {() => (
          <div class={canvasClasses.compactStatus} data-docked="">
            <StatusLine status={status} onOpenLog={props.onOpenLog} />
            <ContextChips status={status} lead={props.context?.()} />
          </div>
        )}
      </Show>
      <Show when={floating.pipe(map((f) => !f))}>
        {() => (
          <DockView
            shell={shell}
            dock={right}
            stripEnd={<Show when={compact}>{() => <RailMenu shell={shell} />}</Show>}
          />
        )}
      </Show>
      <DockView shell={shell} dock={shell.dock('bottom')} />
      {props.commandMenu}
      <BusyOverlay status={status} />
    </div>
  ) as HTMLElement;
  if (props.onDropFiles) {
    root.appendChild(
      (<FileDropOverlay target={root} onDrop={props.onDropFiles} formats={props.dropFormats} />) as Node,
    );
  }
  // Escape closes the open right panel, after anything nearer has had it: handlers on the document (an overlay card, an
  // open list) run first and mark it handled. Not while typing in a field, where Escape belongs to the field.
  const onKeyDown = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented || shell.overlay.value || isEditable(e.target)) return;
    // A pinned panel is part of the workspace, not something passing over the canvas: Escape leaves it.
    if (right.state.value.pinned) return;
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
