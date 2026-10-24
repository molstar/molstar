/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, type ShellModel } from '@molstar/prism-ui';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { onActive } from '@molstar/prism-ui/dom/scope';
import { type ObservableState, observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { CommandManager } from '@molstar/prism-ui/core';
import type { ViewerCommandContext } from '../api/managers.js';
import type { ExtensionUIPort } from '../extensions/host.js';
import type { SceneNode, SceneNodeAction } from '../model/scene.js';
import type { SelectionAction } from '../model/selection.js';
import { SceneRowSlot } from '../views/structure/scene-row-slot.js';
import type { LazyViewerExtension, ViewerExtensionUI, ViewerUIContext } from '../spec.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// Extensions' views in the page's shell (viewer-spec.md §3): sections in the viewer's named slots, tabs,
// cards, pages over the app, rows in the scene tree, and canvas toolbar items. A lazy extension's declared views stand in as placeholders: opening one loads the
// extension, whose real views take their places (same ids).

export function extensionUIPort(
  shell: ShellModel,
  commands: CommandManager<ViewerCommandContext>,
  context: ViewerUIContext,
  /** The selection bar's menu rows the extensions add (`ui.selectionActions`) */
  selectionActions?: ObservableState<readonly SelectionAction[]>,
  /** The scene tree's node menu rows the extensions add (`ui.nodeActions`); absent without a scene */
  nodeActions?: ObservableState<readonly SceneNodeAction[]>,
): ExtensionUIPort {
  // Owners whose load failed: their placeholders, put back by the host, wait for Retry instead of loading again.
  const failed = new Set<string>();
  // A placeholder's body: loading it, or why it did not load and Retry. `then` runs once it loaded.
  const placeholder = (owner: string, load: () => Promise<boolean>, then?: () => void) => {
    const state = observableState<'loading' | 'failed'>(failed.has(owner) ? 'failed' : 'loading');
    const start = () => {
      state.set('loading');
      void load().then((ok) => {
        if (ok) {
          failed.delete(owner);
          then?.();
        } else {
          failed.add(owner);
          state.set('failed');
        }
      });
    };
    // When it is shown, not when it is made: a folded section's or a hidden tab's body is rendered, paused (viewer-spec.md §3.1).
    let started = false;
    onActive((active) => {
      if (!active || started || failed.has(owner)) return;
      started = true;
      start();
    });
    const isFailed = derive(state.asReadonly, (s) => s === 'failed');
    return (
      <div style="padding: var(--msp-space-4); display: grid; gap: var(--msp-space-2); justify-items: start">
        <Show when={isFailed} fallback={() => <p style="opacity: 0.7">Loading…</p>}>
          {() => (
            <>
              <p style="opacity: 0.7">This part could not be loaded.</p>
              <Button size="sm" onClick={start}>
                Retry
              </Button>
            </>
          )}
        </Show>
      </div>
    );
  };
  return {
    context,
    register(_owner, ui: ViewerExtensionUI) {
      const undo: (() => void)[] = [];
      if (ui.tabs?.length) {
        const tab = (t: (typeof ui.tabs)[number]) => ({
          id: t.id,
          title: t.title,
          ...(t.heading ? { heading: t.heading } : {}),
          icon: t.icon,
          ...(t.dock ? { dock: t.dock } : {}),
          ...(t.order !== undefined ? { order: t.order } : {}),
          ...(t.place ? { place: t.place } : {}),
          render: t.render,
        });
        const always = ui.tabs.filter((t) => !t.visible);
        if (always.length) undo.push(shell.registerTabs(always.map(tab)));
        // A tab shown while its `visible` holds: registered then, taken off the rail otherwise.
        for (const t of ui.tabs) {
          if (!t.visible) continue;
          let remove: (() => void) | undefined;
          const sub = t.visible.subscribe((on) => {
            if (on && !remove) remove = shell.registerTabs([tab(t)]);
            else if (!on && remove) {
              remove();
              remove = undefined;
            }
          });
          undo.push(() => {
            sub.unsubscribe();
            remove?.();
          });
        }
      }
      if (ui.overlays?.length) undo.push(shell.registerOverlays(ui.overlays.map((o) => ({ ...o }))));
      if (ui.pages?.length) undo.push(shell.registerPages(ui.pages.map((p) => ({ ...p }))));
      if (ui.sections?.length) {
        undo.push(
          shell.registerSections(
            ui.sections.map((s) => ({
              id: s.id,
              slot: s.slot,
              title: s.title,
              ...(s.heading ? { heading: s.heading(shell) } : {}),
              ...(s.icon ? { icon: s.icon } : {}),
              ...(s.order !== undefined ? { order: s.order } : {}),
              ...(s.keywords ? { keywords: s.keywords } : {}),
              ...(s.open ? { open: true } : {}),
              render: s.render,
            })),
          ),
        );
      }
      // Rows in the scene tree: sections of the tree's slot, shown under the nodes they are for (`SceneTree`).
      if (ui.rows?.length) {
        undo.push(
          shell.registerSections(
            ui.rows.map((r) => ({
              id: r.id,
              slot: SceneRowSlot,
              title: r.id,
              when: (c) => !!c.subject && r.when(c.subject as SceneNode),
              ...(r.whenChanged ? { whenChanged: r.whenChanged } : {}),
              render: r.render,
            })),
          ),
        );
      }
      if (ui.selectionActions?.length && selectionActions) {
        const mine: SelectionAction[] = ui.selectionActions.map((a) => ({
          id: a.id,
          label: a.title,
          icon: a.icon ?? 'mouse-pointer-click',
          ...(a.hint ? { hint: a.hint } : {}),
          ...(a.order !== undefined ? { order: a.order } : {}),
          run: () => (a.run ? a.run() : a.command ? void commands.open(a.command) : undefined),
        }));
        selectionActions.set([...selectionActions.value, ...mine]);
        undo.push(() => selectionActions.set(selectionActions.value.filter((a) => !mine.includes(a))));
      }
      if (ui.nodeActions?.length && nodeActions) {
        const mine: SceneNodeAction[] = ui.nodeActions.map((a) => ({
          id: a.id,
          label: a.title,
          kinds: a.kinds,
          run: a.run,
          ...(a.icon ? { icon: a.icon } : {}),
          ...(a.order !== undefined ? { order: a.order } : {}),
          ...(a.when ? { when: a.when } : {}),
        }));
        nodeActions.set([...nodeActions.value, ...mine]);
        undo.push(() => nodeActions.set(nodeActions.value.filter((a) => !mine.includes(a))));
      }
      if (ui.toolbar?.length) {
        undo.push(
          shell.registerToolbar(
            ui.toolbar.map((t) => ({
              id: t.id,
              icon: t.icon,
              title: t.title,
              ...(t.group ? { group: t.group } : {}),
              ...(t.order !== undefined ? { order: t.order } : {}),
              ...(t.place ? { place: t.place } : {}),
              ...(t.overlay ? { overlay: t.overlay } : {}),
              ...(t.visible ? { visible: t.visible } : {}),
              ...(t.command ? { run: () => void commands.open(t.command!) } : {}),
              ...(t.run ? { run: t.run } : {}),
              ...(t.face ? { face: t.face } : {}),
              ...(t.menu ? { menu: t.menu } : {}),
            })),
          ),
        );
      }
      return () => {
        for (const u of undo.reverse()) u();
      };
    },

    registerStub(owner, stub: LazyViewerExtension['stub'], load) {
      const undo: (() => void)[] = [];
      const tabs = stub.parts.filter((p) => p.kind === 'tab');
      if (tabs.length) {
        undo.push(
          shell.registerTabs(
            tabs.map((t) => ({
              id: t.id,
              title: t.title ?? t.id,
              icon: t.icon ?? 'hourglass',
              ...(t.order !== undefined ? { order: t.order } : {}),
              ...(t.place ? { place: t.place } : {}),
              // Opening it loads the extension; its real tab replaces this one, and is shown (it was asked for).
              render: domView(() => placeholder(owner, load, () => shell.show(t.id))),
            })),
          ),
        );
      }
      const sections = stub.parts.filter((p) => p.kind === 'section' && p.slot);
      if (sections.length) {
        undo.push(
          shell.registerSections(
            sections.map((s) => ({
              id: s.id,
              slot: s.slot!,
              title: s.title ?? s.id,
              ...(s.icon ? { icon: s.icon } : {}),
              ...(s.order !== undefined ? { order: s.order } : {}),
              render: domView(() => placeholder(owner, load)),
            })),
          ),
        );
      }
      const toolbar = stub.parts.filter((p) => p.kind === 'toolbar');
      if (toolbar.length) {
        undo.push(
          shell.registerToolbar(
            toolbar.map((t) => ({
              id: t.id,
              icon: 'hourglass' as const,
              title: t.title ?? t.id,
              run: () => void load(),
            })),
          ),
        );
      }
      return () => {
        for (const u of undo.reverse()) u();
      };
    },
  };
}
