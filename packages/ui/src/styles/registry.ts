/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DomHost } from '../dom/host.js';
import { installTransitions } from './transitions.js';

/**
 * The styles of `@molstar/ui` (styles §9). Each CSS file is compiled into a generated module (`x.css.ts`, see
 * `scripts/css.mjs`) that registers its CSS here when imported, so a page gets exactly the styles of the components it
 * uses. Roots adopt the registered styles into their document or shadow root on connect: no stylesheet to link, in any
 * build setup.
 */

/** Cascade layers, lowest first: host styles outside them always win. */
const LayerOrder = '@layer msui.reset, msui.foundations, msui.primitives, msui.components, msui.overrides;';

const entries = new Map<string, string>();
/** What is written when the registry changes: a document's shared sheet, or a root's `<style>`. */
const targets = new Set<() => void>();
/** A document's one constructed sheet, adopted by it and by its shadow roots: a registration is parsed once. */
const sheets = new Map<Document, { sheet: CSSStyleSheet; update: () => void; users: number }>();
const roots = new Map<StyleRoot, { users: number; remove: () => void }>();

// Registrations come in a burst (a view's modules): apply them once.
let queued = false;
function schedule() {
  if (queued) return;
  queued = true;
  DomHost.queueMicrotask(() => {
    queued = false;
    for (const update of targets) update();
  });
}

/** Registers (or replaces) the CSS of a file; adopted roots update. */
export function defineStyles(id: string, css: string): void {
  if (entries.get(id) === css) return;
  entries.set(id, css);
  if (targets.size) schedule();
}

/** All registered CSS, in registration order, after the layer order. */
export function styleText(): string {
  return [LayerOrder, ...entries.values()].join('\n');
}

/** Where styles are adopted: a document, or a shadow root. */
export type StyleRoot = Document | ShadowRoot;

/**
 * Keeps a document's or a shadow root's styles in sync with the registry: sheets registered later (a lazily loaded
 * view's) reach every adopted root. Adopting a root again adds a user; `releaseStyles` takes one away, and the styles go
 * with the last. A constructed stylesheet where supported, one per document shared by all its roots; otherwise a
 * `<style>` element (in the head, or first in the shadow root).
 */
export function adoptStyles(root: StyleRoot = document): void {
  // The closings go with the styles (styles §7.1): any root that takes the UI's styles, by a `UiElement` or by the app,
  // gets folds that animate closed as they animate open.
  installTransitions();
  const adopted = roots.get(root);
  if (adopted) {
    adopted.users++;
    return;
  }
  const doc = isDocument(root) ? root : root.ownerDocument;
  const win = doc.defaultView;
  let remove: () => void;
  if (win && 'adoptedStyleSheets' in root && typeof win.CSSStyleSheet?.prototype.replaceSync === 'function') {
    let shared = sheets.get(doc);
    if (!shared) {
      const sheet = new win.CSSStyleSheet();
      const update = () => sheet.replaceSync(styleText());
      update();
      shared = { sheet, update, users: 0 };
      sheets.set(doc, shared);
      targets.add(update);
    }
    const docSheet = shared;
    docSheet.users++;
    root.adoptedStyleSheets = [...root.adoptedStyleSheets, docSheet.sheet];
    remove = () => {
      root.adoptedStyleSheets = root.adoptedStyleSheets.filter((s) => s !== docSheet.sheet);
      if (--docSheet.users > 0) return;
      sheets.delete(doc);
      targets.delete(docSheet.update);
    };
  } else {
    const style = doc.createElement('style');
    style.setAttribute('data-msui-styles', '');
    if (isDocument(root)) root.head.appendChild(style);
    else root.prepend(style);
    const update = () => {
      style.textContent = styleText();
    };
    update();
    targets.add(update);
    remove = () => {
      style.remove();
      targets.delete(update);
    };
  }
  roots.set(root, { users: 1, remove });
}

/** Gives up one adoption of a root; with the last, its styles go. It can be adopted again. */
export function releaseStyles(root: StyleRoot): void {
  const adopted = roots.get(root);
  if (!adopted || --adopted.users > 0) return;
  adopted.remove();
  roots.delete(root);
}

function isDocument(root: StyleRoot): root is Document {
  return (root as Document).nodeType === 9;
}
