/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable } from 'rxjs';
import { activity, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { eventTarget, isEditable, ownsKey } from '../../dom/events.js';
import { bind, onCleanup } from '../../dom/scope.js';
import type { CanvasInk } from '../../model/canvas-theme.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { BandDissolvePx, BandFadePx, type ShellModel } from '../../shell/shell.js';
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
  RailToggle,
  StatusLine,
} from './canvas-chrome.js';
import { DockView, FloatingDock } from './dock-view.js';
import { EmptyCanvas } from './empty-canvas.js';
import { FileDropOverlay } from './file-drop.js';
import { PageHost } from './page-host.js';

/**
 * The rail's width, and the one margin of the canvas frame (px, `--msp-canvas-edge`): every floating thing keeps it
 * from the canvas edges, the docked rail, the open panel, and the sequence card.
 */
const RailWidth = 44;
const Edge = 8;
/** What a floating panel leaves clear at the bottom for the corner: the toolbar's row (34px) and the edge gap. */
const CornerRoom = 34 + Edge;

export interface ShellViewProps {
  shell: ShellModel;
  status: StatusModel;
  /** The canvas (the Mol* viewport, or a placeholder in stories). */
  canvas: Child;
  /** Items floating over the canvas besides the toolbar and status line (snapshots, hover label, …). */
  overlays?: Child;
  /** The assistant at work (its turn card, or its line when closed): under the status line's tasks. */
  agent?: Child;
  /** Playback at the top left of the canvas (the snapshot and trajectory card). */
  playback?: Child;
  /**
   * The work in hand at the bottom right of the canvas, before the status chips (the selection bar). A render
   * function: the compact layout docks it elsewhere, and each place renders its own.
   */
  context?: () => Child;
  /**
   * The canvas ink: the theme of the chrome over the canvas, from the canvas background (light background: light glass,
   * independent of the UI theme: dark glass over a white canvas would be muddy), and the background color the chrome is
   * tinted from instead of the UI's own background, so cards read as frosted canvas rather than slabs of panel color.
   */
  ink?: ReadonlyState<CanvasInk>;
  /** The command menu while it is open: over the whole shell. */
  commandMenu?: Child;
  /** The top of the canvas stack, above the canvas toolbar items: the interact / select mode buttons. */
  canvasLead?: Child;
  /** Beside the canvas stack, outside its frame: what opens from it (the picking level, next to Select). */
  canvasAside?: Child;
  onOpenLog?: () => void;
  /** Nothing loaded: a while after, the Mol* wordmark shows in the middle of the canvas. */
  empty?: Observable<boolean>;
  /** The splash under the mark once the viewer has been idle on the empty canvas a while (on by default). */
  splash?: boolean;
  /** The viewer is working: the splash waits for it to be idle. */
  busy?: Observable<boolean>;
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
  // Expanded (`view.expand`): the canvas fills the viewer; the rail, the docks and the band are left out (not rendered),
  // and the insets they take go to the frame's margin.
  const expanded = shell.expanded;
  const shown = expanded.pipe(map((e) => !e));
  // The right dock floats over the canvas unless it is docked, or the shell is compact (then it is a bottom sheet).
  const floating = combineLatest([right.state, shell.mode]).pipe(
    map(([s, mode]) => s.floating && mode !== 'compact'),
    distinctUntilChanged(),
  );
  // The top dock (the sequence band) floats in the compact layout too: a band of the canvas, not a panel above it.
  const floatingTop = top.state.pipe(
    map((s) => s.floating),
    distinctUntilChanged(),
  );
  // A top tab with no frame is a band of the canvas area (the sequence), over the canvas' top: the canvas' own color
  // behind its lines, fading out below them (BandFadePx); the canvas' chrome goes below the lines.
  // What floats at the top (the sequence card, or the text of a band) takes this much: the snapshots and the canvas stack
  // go below it (a band's faded bottom is not in the way).
  const topInset = combineLatest([top.state, top.active, floatingTop, expanded]).pipe(
    map(([s, tab, f, e]) =>
      f && tab && !e ? (tab.frame === 'none' ? s.size - BandFadePx + Edge : Edge + s.size + Edge) : Edge,
    ),
    distinctUntilChanged(),
    map((px) => `${px}px`),
  );
  // A compact rail's column hidden: the canvas takes its width.
  const railHidden = combineLatest([floating, shell.railShown, expanded]).pipe(
    map(([f, shown, e]) => (f && !shown) || e),
    distinctUntilChanged(),
  );
  const railPx = railHidden.pipe(map((h) => (h ? 0 : RailWidth)));
  // A compact rail lies over the canvas' edge (it does not take the canvas' width, so opening it moves nothing): the
  // canvas runs to the edge, and the corner keeps clear of the rail's column (where its toggle sits) whether it shows.
  const compactRail = combineLatest([floating, shell.compactRail, expanded]).pipe(
    map(([f, c, e]) => f && c && !e),
    distinctUntilChanged(),
  );
  const railTake = combineLatest([railPx, compactRail]).pipe(
    map(([px, c]) => (c ? 0 : px)),
    distinctUntilChanged(),
  );
  // The button that shows and hides it, in the canvas's corner: only where the rail floats, and the rail is compact.
  const railToggle = combineLatest([floating, shell.compactRail, expanded]).pipe(
    map(([f, c, e]) => f && c && !e),
    distinctUntilChanged(),
  );
  // Pinned: the canvas ends where the rail, or the open panel, begins (they sit beside it, not over it).
  const pinned = combineLatest([right.state, floating, expanded]).pipe(
    map(([s, f, e]) => f && s.pinned && !e),
    distinctUntilChanged(),
  );
  // An overlay card open in the panel's place counts as the panel.
  const pinnedWidth = combineLatest([right.state, right.active, shell.overlay, pinned, railPx]).pipe(
    map(([s, active, overlay, p, rail]) => (p ? rail + (active || overlay ? s.size : 0) : 0)),
    distinctUntilChanged(),
  );
  // The canvas and the band end where the fixed part of the layout begins: a pinned panel and the rail, or the rail alone
  // while the panel floats (it lies over both).
  const canvasEnd = combineLatest([pinnedWidth, floating, railTake]).pipe(
    map(([pinnedPx, f, rail]) => (pinnedPx > 0 ? pinnedPx : f ? rail : 0)),
    distinctUntilChanged(),
  );
  // What floats on the right takes this much of the canvas: the bottom bar keeps clear of it.
  const rightInset = combineLatest([right.state, right.active, floating, pinned, railPx, expanded]).pipe(
    // The docked rail is the canvas's right edge, and the open panel is too: what sits beside either keeps the frame
    // margin. Pinned, both are outside the canvas: only the margin is left.
    map(([s, active, f, p, rail, e]) => (p || e ? Edge : rail + Edge + (f && active ? s.size + Edge : 0))),
    distinctUntilChanged(),
    map((px) => `${px}px`),
  );
  const cornerInset = combineLatest([pinned, compact, railTake, compactRail, expanded]).pipe(
    map(([p, c, rail, overlay, e]) => (p || c || e ? Edge : (overlay ? RailWidth : rail) + Edge)),
    distinctUntilChanged(),
    map((px) => `${px}px`),
  );
  // An open floating right panel: overlays then open in its place instead of stacking glass on glass.
  const rightPanelOpen = combineLatest([right.active, floating, expanded]).pipe(
    map(([active, f, e]) => (f && active && !e ? 'open' : undefined)),
    distinctUntilChanged(),
  );
  const rightPanelWidth = right.state.pipe(map((s) => `${s.size}px`));
  const root = (
    <div class={classes.shell} data-mode={shell.mode} data-expanded={expanded.pipe(map((e) => (e ? '' : undefined)))}>
      <Show when={combineLatest([floatingTop, expanded]).pipe(map(([f, e]) => !f && !e))}>
        {() => <DockView shell={shell} dock={top} />}
      </Show>
      <Show when={shown}>{() => <DockView shell={shell} dock={shell.dock('left')} />}</Show>
      <main
        class={classes.viewport}
        // The canvas and what floats over it (the sequence band, cards): the hover label follows the pointer here.
        data-msp-viewport=""
        style={{
          '--msp-canvas-edge': `${Edge}px`,
          '--msp-canvas-right-inset': rightInset,
          // Where the bottom-right corner (the toolbar and the selection) ends: it does not move with a floating panel,
          // which stops above it (`--msp-corner-room`); pinned, the canvas itself ends at the panel.
          '--msp-canvas-corner-inset': cornerInset,
          '--msp-corner-room': `${CornerRoom}px`,
          '--msp-canvas-top-inset': topInset,
          // Where the canvas and a band end, px from the right: the rail (and a pinned panel). Anchored to the viewport's
          // right edge.
          '--msp-canvas-end': canvasEnd.pipe(map((px) => `${px}px`)),
          '--msp-band-fade': `${BandFadePx}px`,
          '--msp-band-dissolve': `${BandDissolvePx}px`,
          // Clear of the rail only (what sits beside it whatever the panel does: the compact stack).
          '--msp-canvas-rail-inset': railPx.pipe(map((px) => `${px + Edge}px`)),
          '--msp-right-panel-width': rightPanelWidth,
          // Where cards beside the rail end (the floating panel's right edge).
          '--msp-canvas-panel-right': railPx.pipe(map((px) => `${px + Edge}px`)),
          // How much of the right the pinned rail and panel take from the canvas (0 when floating over it).
          '--msp-canvas-pinned': pinnedWidth.pipe(map((px) => `${px}px`)),
        }}
      >
        <div class={classes.canvas} data-msp-canvas="">
          {props.canvas}
        </div>
        {/* Chrome that is part of the view: themed by the canvas background, tinted from its color. */}
        <div
          class={`${classes.layer} ${classes.canvasLayer}`}
          data-msp-theme={props.ink ? props.ink.pipe(map((i) => i.theme)) : 'light'}
          style={{
            // Every derived shade (glass, controls, hover, muted text) follows from the base color.
            '--msp-color-bg': props.ink ? props.ink.pipe(map((i) => i.color)) : undefined,
          }}
        >
          <Show when={combineLatest([floatingTop, expanded]).pipe(map(([f, e]) => f && !e))}>
            {() => <FloatingDock shell={shell} dock={top} />}
          </Show>
          {props.empty ? (
            <EmptyCanvas
              empty={props.empty}
              splash={props.splash ?? true}
              {...(props.busy ? { busy: props.busy } : {})}
            />
          ) : null}
          <div class={canvasClasses.lead}>{props.playback}</div>
          {props.overlays}
          {/* The compact rail's button: in the rail's column, at the bottom, where the rail slides in around it. */}
          <Show when={railToggle}>{() => <RailToggle shell={shell} />}</Show>
          {/* The bottom of the canvas: the status on the left, and in the corner on the right the work in hand (the
              selection) to the left of the canvas toolbar on its row. Compact: the status and the selection are docked under the
              canvas; the toolbar stays in the corner. */}
          <div class={canvasClasses.bottom}>
            <div class={canvasClasses.bottomGrid}>
              <Show when={compact.pipe(map((c) => !c))}>
                {() => <StatusLine status={status} onOpenLog={props.onOpenLog} after={props.agent} />}
              </Show>
              <div class={canvasClasses.corner}>
                <CanvasStack shell={shell} lead={props.canvasLead} aside={props.canvasAside} />
                <Show when={compact.pipe(map((c) => !c))}>
                  {() => <ContextChips status={status} lead={props.context?.()} />}
                </Show>
              </div>
            </div>
          </div>
        </div>
        {/* The application's own chrome over the canvas (the rail, the panel, the cards): the UI theme. */}
        <div
          class={classes.layer}
          data-right-panel={rightPanelOpen}
          data-pinned={pinned.pipe(map((p) => (p ? '' : undefined)))}
          data-overlay={shell.overlay.pipe(map((o) => (o ? 'open' : undefined)))}
        >
          <OverlayHost shell={shell} />
          {/* Over the canvas, beside the rail and the open panel (they stay, and the page makes room for them). */}
          <PageHost shell={shell} />
          <Show when={shown}>
            {() => (
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
                {() => {
                  // A compact rail's hidden column is inert and paused (`activity`); what hides it is bound out here, since
                  // what the activity renders pauses with it.
                  const dock = activity(shell.railShown, () => (
                    <FloatingDock
                      shell={shell}
                      dock={right}
                      railItems={(rail) => (
                        <CanvasToolbar
                          shell={shell}
                          part="rest"
                          hidden={rail.hidden}
                          railMenu={(item, anchor) => rail.toggleMenu(anchor, item)}
                        />
                      )}
                      railFoot={(rail) => (
                        <CanvasToolbar
                          shell={shell}
                          part="lead"
                          hidden={rail.hidden}
                          railMenu={(item, anchor) => rail.toggleMenu(anchor, item)}
                        />
                      )}
                    />
                  )) as HTMLElement;
                  bind(shell.railShown, (shown) => dock.toggleAttribute('data-rail-hidden', !shown));
                  bind(shell.compactRail, (c) => dock.toggleAttribute('data-rail-overlay', c));
                  return dock;
                }}
              </Show>
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
      <Show when={combineLatest([floating, expanded]).pipe(map(([f, e]) => !f && !e))}>
        {() => (
          <DockView
            shell={shell}
            dock={right}
            stripEnd={<Show when={compact}>{() => <RailMenu shell={shell} />}</Show>}
          />
        )}
      </Show>
      <Show when={shown}>{() => <DockView shell={shell} dock={shell.dock('bottom')} />}</Show>
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
    if (e.key !== 'Escape' || e.defaultPrevented || shell.overlay.value || isEditable(eventTarget(e))) return;
    if (!ownsKey(root, e)) return;
    // Expanded: Escape brings the chrome back (`view.expand` toggled off).
    if (expanded.value) {
      shell.setExpanded(false);
      e.preventDefault();
      return;
    }
    // A pinned panel is part of the workspace, not something passing over the canvas: Escape leaves it.
    if (right.state.value.pinned) return;
    if (!right.state.value.collapsed && right.active.value !== undefined) {
      shell.toggle('right');
      e.preventDefault();
    }
  };
  document.addEventListener('keydown', onKeyDown);
  onCleanup(() => document.removeEventListener('keydown', onKeyDown));
  // Tooltips for the icon buttons of the rail and the canvas stack.
  onCleanup(installTooltips(root));
  // The shell decides compact or desktop by its own width, not the window's.
  const observer = new ResizeObserver((entries) => shell.setWidth(entries[0].contentRect.width));
  observer.observe(root);
  onCleanup(() => observer.disconnect());
  return root;
}
