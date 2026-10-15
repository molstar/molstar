/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { MessageKind } from '@molstar/prism-ui/core';

// The plugin's toasts (`plugin.managers.toast`: what behaviors and extensions show with `PluginCommands.Toast.Show`)
// as the viewer's status messages (parity.md's "toasts, task progress" row; journal.md §3): each new toast, or one
// shown again under its key, is posted once. The status line fades it as its own messages; the Log keeps it.

export interface ToastPort {
  post(text: string, kind?: MessageKind): void;
}

/** Posts the plugin's toasts to `notifications`; returns the unsubscribe. */
export function watchToasts(plugin: PluginContext, notifications: ToastPort): () => void {
  const toasts = plugin.managers.toast;
  // The entries already posted, by id (a toast shown again under its key gets a new id).
  const posted = new Set<number>();
  const sub = toasts.events.changed.subscribe(() => {
    const live = new Set<number>();
    toasts.state.entries.forEach((e) => {
      if (!e) return;
      live.add(e.id);
      if (posted.has(e.id)) return;
      posted.add(e.id);
      const text = toastText(e.title, e.message);
      if (text) notifications.post(text, 'info');
    });
    // Forget the ones that went: the set holds only what is shown.
    for (const id of posted) if (!live.has(id)) posted.delete(id);
  });
  return () => sub.unsubscribe();
}

/** A toast in words: its title and its message as text (an HTML string's text; a component's title only). */
export function toastText(title: string, message: string | Function): string {
  const body = typeof message === 'string' ? htmlText(message) : '';
  if (!body) return title.trim();
  if (!title.trim()) return body;
  return `${title.trim()}: ${body}`;
}

function htmlText(html: string): string {
  return html
    .replace(/<br\s*\/?>/gi, ' ')
    .replace(/<[^>]*>/g, '')
    .replace(/&nbsp;/g, ' ')
    .replace(/&lt;/g, '<')
    .replace(/&gt;/g, '>')
    .replace(/&quot;/g, '"')
    .replace(/&#39;/g, "'")
    .replace(/&amp;/g, '&')
    .replace(/\s+/g, ' ')
    .trim();
}
