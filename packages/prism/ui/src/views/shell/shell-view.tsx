/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { activity, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { eventTarget, isEditable, ownsKey } from '../../dom/events.js';
import { type MaybeReactive, onCleanup, type WatchSource, watch } from '../../dom/scope.js';
import { observeResize, on } from '../../dom/scoped.js';
import { setAttr, setStyle } from '../../dom/writes.js';
import type { CanvasInk } from '../../model/canvas-theme.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
import { BandDissolvePx, BandFadePx, type CanvasMenuItem, type ShellModel } from '../../shell/shell.js';
import type { NotificationsManager } from '../../managers/notifications.js';
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
import { DockView, dockFacets, FloatingDock } from './dock-view.js';
import { EmptyCanvas } from './empty-canvas.js';
import type { SplashKind } from './prism-splash.js';
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
  status: NotificationsManager;
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
  /**
   * The rows of `canvasLead` for the collapsed toolbar's menu (the mode choices, a submenu for their options): the
   * stack lists them above the toolbar items (shell.md §2). Without it the menu lists the items only.
   */
  canvasLeadMenu?: MaybeReactive<readonly CanvasMenuItem[]>;
  /** Beside the canvas stack, outside its frame: what opens from it (the picking level, next to Select). */
  canvasAside?: Child;
  onOpenLog?: () => void;
  /** Nothing loaded: a while after, the Mol* wordmark shows in the middle of the canvas. */
  empty?: MaybeReactive<boolean>;
  /**
   * The splash under the mark once the viewer has been idle on the empty canvas a while (on by default), and which:
   * the light's particles, the light as a wave, or either at random each time it shows (the default).
   */
  splash?: boolean;
  splashField?: SplashKind | 'random';
  /**
   * The viewer is working now (not debounced: a host's operation, a chunk loading, its engine updating): the splash waits
   * for it to be idle. A debounced flag comes later than the splash's wait and lets it show during a load.
   */
  busy?: MaybeReactive<boolean>;
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
  const compact = derive(shell.mode, (m) => m === 'compact');
  const top = shell.dock('top');
  // Expanded (`view.expand`): the canvas fills the viewer; the rail, the docks and the band are left out (not rendered),
  // and the insets they take go to the frame's margin.
  const expanded = shell.expanded;
  const shown = derive(expanded, (e) => !e);
  // The right dock floats over the canvas unless it is docked, or the shell is compact (then it is a bottom sheet).
  // What the right dock's state decides here, through its facets (a change of its state wakes each once).
  const rightDock = dockFacets(shell, right);
  const floating = derive([rightDock.floating, shell.mode], (f, mode) => f && mode !== 'compact');
  // The top dock (the sequence band) floats in the compact layout too: a band of the canvas, not a panel above it.
  const floatingTop = dockFacets(shell, top).floating;
  // A top tab with no frame is a band of the canvas area (the sequence), over the canvas' top: the canvas' own color
  // behind its lines, fading out below them (BandFadePx); the canvas' chrome goes below the lines.
  // What floats at the top (the sequence card, or the text of a band) takes this much: the snapshots and the canvas stack
  // go below it (a band's faded bottom is not in the way).
  const topInset = derive(
    [dockFacets(shell, top).size, top.active, floatingTop, expanded],
    (size, tab, f, e) =>
      `${f && tab && !e ? (tab.frame === 'none' ? size - BandFadePx + Edge : Edge + size + Edge) : Edge}px`,
  );
  // A compact rail's column hidden: the canvas takes its width.
  const railPx = derive([floating, shell.railShown, expanded], (f, shown, e) => ((f && !shown) || e ? 0 : RailWidth));
  // A compact rail lies over the canvas' edge (it does not take the canvas' width, so opening it moves nothing): the
  // canvas runs to the edge, and the corner keeps clear of the rail's column (where its toggle sits) whether it shows.
  const compactRail = derive([floating, shell.compactRail, expanded], (f, c, e) => f && c && !e);
  const railTake = derive([railPx, compactRail], (px, c) => (c ? 0 : px));
  // The button that shows and hides it, in the canvas's corner: only where the rail floats, and the rail is compact.
  const railToggle = compactRail;
  // Pinned: the canvas ends where the rail, or the open panel, begins (they sit beside it, not over it).
  const pinned = derive([rightDock.pinned, floating, expanded], (p, f, e) => f && p && !e);
  // An overlay card open in the panel's place counts as the panel.
  const pinnedWidth = derive(
    [rightDock.size, right.active, shell.overlay, pinned, railPx],
    (size, active, overlay, p, rail) => (p ? rail + (active || overlay ? size : 0) : 0),
  );
  // The canvas and the band end where the fixed part of the layout begins: a pinned panel and the rail, or the rail alone
  // while the panel floats (it lies over both).
  const canvasEnd = derive([pinnedWidth, floating, railTake], (pinnedPx, f, rail) =>
    pinnedPx > 0 ? pinnedPx : f ? rail : 0,
  );
  // What floats on the right takes this much of the canvas: the bottom bar and the hover label keep clear of it.
  const rightInset = derive(
    [rightDock.size, right.active, shell.overlay, floating, pinned, railPx, expanded, compact],
    // The docked rail is the canvas's right edge, and the open panel is too, or an overlay card open in its place: what
    // sits beside either keeps the frame margin. Pinned, both are outside the canvas: only the margin is left. Compact,
    // the rail is the bottom bar.
    (size, active, overlay, f, p, rail, e, c) =>
      `${p || e ? Edge : (c ? 0 : rail) + Edge + (f && (active || !!overlay) ? size + Edge : 0)}px`,
  );
  const cornerInset = derive(
    [pinned, compact, railTake, compactRail, expanded],
    (p, c, rail, overlay, e) => `${p || c || e ? Edge : (overlay ? RailWidth : rail) + Edge}px`,
  );
  // An open floating right panel: overlays then open in its place instead of stacking glass on glass.
  const rightPanelOpen = derive([right.active, floating, expanded], (active, f, e) =>
    f && active && !e ? 'open' : undefined,
  );
  const rightPanelWidth = rightDock.sizePx;
  const railInset = derive(railPx, (px) => `${px + Edge}px`);
  // The canvas: quiet (transparent, its axes and all) from the start and while the empty canvas's mark and splash are
  // up, coming in gradually over its own color when something is in it (shell.module.css `.canvas[data-quiet]`); a
  // load the mark never showed for finds it shown. Keys scoped to the canvas work while the focus is in it
  // (`data-msp-scope`, commands.md §6).
  const canvasBox = (
    <div class={classes.canvas} data-msp-canvas="" data-msp-scope="canvas" data-quiet={props.empty ? '' : undefined}>
      {props.canvas}
    </div>
  ) as HTMLElement;
  const canvasShown = (shown: boolean) => canvasBox.toggleAttribute('data-quiet', !shown);
  let viewport!: HTMLElement;
  let canvasLayer!: HTMLElement;
  let chromeLayer!: HTMLElement;
  const root = (
    <div class={classes.shell}>
      <Show when={derive([floatingTop, expanded], (f, e) => !f && !e)}>
        {() => <DockView shell={shell} dock={top} />}
      </Show>
      <Show when={shown}>{() => <DockView shell={shell} dock={shell.dock('left')} />}</Show>
      <main
        class={classes.viewport}
        // The canvas and what floats over it (the sequence band, cards): the hover label follows the pointer here.
        data-msp-viewport=""
        // What never changes; what follows the layout is written by the shell's watches below.
        style={{
          '--msp-canvas-edge': `${Edge}px`,
          '--msp-corner-room': `${CornerRoom}px`,
          '--msp-band-fade': `${BandFadePx}px`,
          '--msp-band-dissolve': `${BandDissolvePx}px`,
        }}
        ref={(el) => (viewport = el)}
      >
        {canvasBox}
        {/* Chrome that is part of the view: themed by the canvas background, tinted from its color. */}
        <div
          class={`${classes.layer} ${classes.canvasLayer}`}
          data-msp-theme={props.ink ? undefined : 'light'}
          ref={(el) => (canvasLayer = el)}
        >
          <Show when={derive([floatingTop, expanded], (f, e) => f && !e)}>
            {() => <FloatingDock shell={shell} dock={top} />}
          </Show>
          {props.empty ? (
            <EmptyCanvas
              empty={props.empty}
              splash={props.splash ?? true}
              field={props.splashField ?? 'random'}
              onCanvasShown={canvasShown}
              {...(props.busy ? { busy: props.busy } : {})}
            />
          ) : null}
          <div class={canvasClasses.lead}>{props.playback}</div>
          {props.overlays}
          {/* The compact rail's button: in the rail's column, at the bottom, where the rail slides in around it. */}
          <Show when={railToggle}>{() => <RailToggle shell={shell} />}</Show>
          {/* The bottom of the canvas: the status on the left, and in the corner on the right the work in hand (the
              selection) to the left of the canvas toolbar on its row. Compact too: on the canvas, so a message or a selection
              coming and going moves nothing. */}
          {/* It acts on the scene as the canvas does: a press on it leaves an open card open (`OverlayCard`). */}
          <div class={canvasClasses.bottom} data-keeps-card="">
            <div class={canvasClasses.bottomGrid}>
              <StatusLine status={status} onOpenLog={props.onOpenLog} after={props.agent} />
              <div class={canvasClasses.corner}>
                <CanvasStack
                  shell={shell}
                  lead={props.canvasLead}
                  aside={props.canvasAside}
                  {...(props.canvasLeadMenu ? { leadMenu: props.canvasLeadMenu } : {})}
                  quiet={props.empty}
                />
                <ContextChips status={status} lead={props.context?.()} />
              </div>
            </div>
          </div>
        </div>
        {/* The application's own chrome over the canvas (the rail, the panel, the cards): the UI theme. */}
        <div class={classes.layer} ref={(el) => (chromeLayer = el)}>
          <OverlayHost shell={shell} />
          {/* Over the canvas, beside the rail and the open panel (they stay, and the page makes room for them). */}
          <PageHost shell={shell} />
          <Show when={shown}>
            {() => (
              <Show
                when={floating}
                fallback={() => (
                  // Docked beside the canvas (not compact: there the rail's actions are under More).
                  <Show when={derive(compact, (c) => !c)}>
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
                  watch([shell.railShown, shell.compactRail], () => {
                    setAttr(dock, 'data-rail-hidden', !shell.railShown.value);
                    setAttr(dock, 'data-rail-overlay', shell.compactRail.value);
                  });
                  return dock;
                }}
              </Show>
            )}
          </Show>
        </div>
      </main>
      <Show when={derive([floating, expanded], (f, e) => !f && !e)}>
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
  // The frame as a function of the layout (each source a gate over the dock and shell states, above): what floats takes
  // its insets from these; one watch for the viewport's variables, others for the attributes of the layers.
  watch([shell.mode, expanded], () => {
    setAttr(root, 'data-mode', shell.mode.value);
    setAttr(root, 'data-expanded', expanded.value);
  });
  const insets: WatchSource[] = [rightInset, cornerInset, topInset, canvasEnd, railInset, rightPanelWidth, pinnedWidth];
  if (props.ink) insets.push(props.ink);
  watch(insets, () => {
    setStyle(viewport, '--msp-canvas-right-inset', rightInset.value);
    // Where the bottom-right corner (the toolbar and the selection) ends: it does not move with a floating panel, which
    // stops above it (`--msp-corner-room`); pinned, the canvas itself ends at the panel.
    setStyle(viewport, '--msp-canvas-corner-inset', cornerInset.value);
    setStyle(viewport, '--msp-canvas-top-inset', topInset.value);
    // Where the canvas and a band end, px from the right: the rail (and a pinned panel). Anchored to the viewport's
    // right edge.
    setStyle(viewport, '--msp-canvas-end', `${canvasEnd.value}px`);
    // Clear of the rail only (what sits beside it whatever the panel does: the compact stack); where cards beside the
    // rail end (the floating panel's right edge).
    setStyle(viewport, '--msp-canvas-rail-inset', railInset.value);
    setStyle(viewport, '--msp-canvas-panel-right', railInset.value);
    setStyle(viewport, '--msp-right-panel-width', rightPanelWidth.value);
    // How much of the right the pinned rail and panel take from the canvas (0 when floating over it).
    setStyle(viewport, '--msp-canvas-pinned', `${pinnedWidth.value}px`);
    // The canvas's own color under it: a canvas fading in (`reveal`) shows the structure coming, not the shell.
    setStyle(viewport, 'background', props.ink?.value.color);
  });
  if (props.ink) {
    const ink = props.ink;
    watch([ink], () => {
      setAttr(canvasLayer, 'data-msp-theme', ink.value.theme);
      // Every derived shade (glass, controls, hover, muted text) follows from the base color.
      setStyle(canvasLayer, '--msp-color-bg', ink.value.color);
    });
  }
  watch([rightPanelOpen, pinned, shell.overlay], () => {
    setAttr(chromeLayer, 'data-right-panel', rightPanelOpen.value);
    setAttr(chromeLayer, 'data-pinned', pinned.value);
    setAttr(chromeLayer, 'data-overlay', shell.overlay.value ? 'open' : null);
  });
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
  on(document, 'keydown', onKeyDown);
  // Tooltips for the icon buttons of the rail and the canvas stack.
  onCleanup(installTooltips(root));
  // The shell decides compact or desktop by its own width, not the window's.
  observeResize(root, (entries) => shell.setWidth(entries[0].contentRect.width));
  return root;
}
