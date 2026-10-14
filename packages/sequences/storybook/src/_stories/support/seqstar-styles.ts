/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stylesheets of the Seq* stories. Prism's tokens and reset are prism-ui's own (`@molstar/prism-ui/styles/base`),
// adopted into the document as Prism's elements do; the host Storybook adopts them for Prism's stories already, and the
// adoption is idempotent. The page defaults (`page.css`) are Seq*'s: the host shows other packages' stories in the same
// page (iframe), so they never go in globally. A story page takes them while it is mounted; the last to leave removes
// them.

import '@molstar/prism-ui/styles/base';
import { adoptStyles } from '@molstar/prism-ui/styles/registry';
import { useInsertionEffect } from 'react';
import pageCss from '../ui/page.css?inline';

const SHEETS = { page: pageCss } as const;
type Sheet = 'tokens' | keyof typeof SHEETS;

const held = new Map<keyof typeof SHEETS, { style: HTMLStyleElement; users: number }>();

function acquire(sheet: Sheet): () => void {
  if (sheet === 'tokens') {
    adoptStyles(document);
    return () => {};
  }
  let h = held.get(sheet);
  if (!h) {
    const style = document.createElement('style');
    style.setAttribute('data-seq-stories', sheet);
    style.textContent = SHEETS[sheet];
    document.head.appendChild(style);
    h = { style, users: 0 };
    held.set(sheet, h);
  }
  const entry = h;
  entry.users++;
  return () => {
    if (--entry.users > 0) return;
    entry.style.remove();
    held.delete(sheet);
  };
}

/** The stylesheets a component needs while it is mounted (before it first lays out). */
export function useSeqStyles(...sheets: Sheet[]) {
  const key = sheets.join();
  useInsertionEffect(() => {
    const releases = sheets.map(acquire);
    return () => {
      for (const release of releases) release();
    };
  }, [key]);
}
