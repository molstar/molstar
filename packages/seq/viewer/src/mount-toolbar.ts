/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Toolbar, type ToolbarItem } from '@molstar/ui/components/toolbar/toolbar';
import { installTooltips } from '@molstar/ui/components/tooltip/tooltip';
import { onCleanup, render } from '@molstar/ui/dom';
import { DefaultToolbar, toolbarButtons } from './toolbar.js';
import type { SeqViewer } from './viewer.js';

// The mounted toolbar (_spec/README.md "Input", Toolbar): `@molstar/ui`'s toolbar over `toolbarButtons`, in a block of its
// own above the canvas. Everything it adds is taken away by the function it returns.

/** Mounts the bar for the commands named (`ToolbarMount`), above whatever is in `target`. */
export function mountToolbar(
  viewer: SeqViewer,
  target: HTMLElement,
  commands: readonly string[] = DefaultToolbar,
): () => void {
  const doc = target.ownerDocument;
  const host = doc.createElement('div');
  host.style.position = 'relative';
  host.className = 'seqx-toolbar';
  // the toolbar's own disabled colour (text-subtle) is a few percent from its enabled one (text-muted): too near to
  // tell a disabled button from an enabled one, in either theme; opacity separates them in both
  const style = doc.createElement('style');
  style.textContent = '.seqx-toolbar [aria-disabled="true"] { opacity: 0.4; }';
  host.append(style);
  // before the canvas the viewer mounted into the target
  target.prepend(host);
  const items: ToolbarItem[] = toolbarButtons(viewer, commands).map((b) => ({
    id: b.command,
    icon: b.icon as ToolbarItem['icon'],
    title: b.title,
    disabled: b.disabled,
    run: () => void b.run(),
  }));
  // the tooltips are the render scope's: they go with it
  const dispose = render(() => {
    onCleanup(installTooltips(host));
    return Toolbar({ groups: [{ id: 'view', label: 'View', items }], 'aria-label': 'Sequence view' });
  }, host);
  return () => {
    dispose();
    host.remove();
  };
}
