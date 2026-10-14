/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stylesheets of the Seq* stories, in the document only while something of Seq* uses them. The host Storybook shows
// other packages' stories in the same page (iframe), so these never go in globally: a story page takes the tokens and the
// page defaults, a docs example the tokens; the last to leave removes them.

import { useInsertionEffect } from 'react';
import pageCss from '../ui/page.css?inline';
import tokensCss from '../ui/tokens.css?inline';

const SHEETS = { tokens: tokensCss, page: pageCss } as const;
type Sheet = keyof typeof SHEETS;

const held = new Map<Sheet, { style: HTMLStyleElement; users: number }>();

function acquire(sheet: Sheet): () => void {
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
