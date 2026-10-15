/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, type ShellModel } from '@molstar/prism-ui';
import type { CommandManager } from '@molstar/prism-ui/core';
import type { ViewerCommandContext } from '../api/managers.js';
import type { ExtensionUIPort } from '../extensions/host.js';
import type { LazyViewerExtension, ViewerExtensionUI, ViewerUIContext } from '../spec.js';

// Extensions' views in the page's shell (viewer-spec.md §3, README §3.1): sections in the viewer's named slots, tabs,
// cards, and canvas toolbar items. A lazy extension's declared views stand in as placeholders: opening one loads the
// extension, whose real views take their places (same ids).

export function extensionUIPort(
  shell: ShellModel,
  commands: CommandManager<ViewerCommandContext>,
  context: ViewerUIContext,
): ExtensionUIPort {
  return {
    context,
    register(_owner, ui: ViewerExtensionUI) {
      const undo: (() => void)[] = [];
      if (ui.tabs?.length) {
        undo.push(
          shell.registerTabs(
            ui.tabs.map((t) => ({
              id: t.id,
              title: t.title,
              icon: t.icon,
              ...(t.dock ? { dock: t.dock } : {}),
              ...(t.order !== undefined ? { order: t.order } : {}),
              render: t.render,
            })),
          ),
        );
      }
      if (ui.overlays?.length) undo.push(shell.registerOverlays(ui.overlays.map((o) => ({ ...o }))));
      if (ui.sections?.length) {
        undo.push(
          shell.registerSections(
            ui.sections.map((s) => ({
              id: s.id,
              slot: s.slot,
              title: s.title,
              ...(s.icon ? { icon: s.icon } : {}),
              ...(s.order !== undefined ? { order: s.order } : {}),
              ...(s.keywords ? { keywords: s.keywords } : {}),
              render: s.render,
            })),
          ),
        );
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
            })),
          ),
        );
      }
      return () => {
        for (const u of undo.reverse()) u();
      };
    },

    registerStub(_owner, stub: LazyViewerExtension['stub'], load) {
      const undo: (() => void)[] = [];
      const loading = () => <p style="padding: var(--msp-space-4); opacity: 0.7">Loading…</p>;
      const tabs = stub.parts.filter((p) => p.kind === 'tab');
      if (tabs.length) {
        undo.push(
          shell.registerTabs(
            tabs.map((t) => ({
              id: t.id,
              title: t.title ?? t.id,
              icon: t.icon ?? 'hourglass',
              ...(t.order !== undefined ? { order: t.order } : {}),
              // Opening it loads the extension; its real tab replaces this one, and is shown (it was asked for).
              render: domView(() => {
                void load().then((ok) => {
                  if (ok) shell.show(t.id);
                });
                return loading();
              }),
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
              render: domView(() => {
                void load();
                return loading();
              }),
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
