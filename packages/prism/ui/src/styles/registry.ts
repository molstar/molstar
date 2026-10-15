/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The styles of prism-ui. Each CSS file is compiled into a generated module (`x.css.ts`, see `scripts/css.mjs`) that
 * registers its CSS here when imported, so a page gets exactly the styles of the components it uses. Elements adopt
 * the registered styles into their document on connect: no stylesheet to link, in any build setup.
 */

/** Cascade layers, lowest first: host styles outside them always win. */
const LayerOrder = '@layer msp.reset, msp.theme, msp.components, msp.overrides;';

const entries = new Map<string, string>();
const roots = new Map<StyleRoot, () => void>();

/** Registers (or replaces) the CSS of a file; adopted documents update. */
export function defineStyles(id: string, css: string): void {
  if (entries.get(id) === css) return;
  entries.set(id, css);
  for (const update of roots.values()) update();
}

/** All registered CSS, in registration order, after the layer order. */
export function styleText(): string {
  return [LayerOrder, ...entries.values()].join('\n');
}

/** Where styles are adopted: a document, or a shadow root (a `<prism-viewer>`'s). */
export type StyleRoot = Document | ShadowRoot;

/**
 * Keeps a document's or a shadow root's styles in sync with the registry (idempotent): sheets registered later (a lazy
 * extension's) reach every adopted root. Uses a constructed stylesheet where supported, otherwise a `<style>` element
 * (in the head, or first in the shadow root).
 */
export function adoptStyles(root: StyleRoot = document): void {
  if (roots.has(root)) return;
  const doc = isDocument(root) ? root : root.ownerDocument;
  const win = doc.defaultView;
  let update: () => void;
  if (win && 'adoptedStyleSheets' in root && typeof win.CSSStyleSheet?.prototype.replaceSync === 'function') {
    const sheet = new win.CSSStyleSheet();
    root.adoptedStyleSheets = [...root.adoptedStyleSheets, sheet];
    update = () => sheet.replaceSync(styleText());
  } else {
    const style = doc.createElement('style');
    style.setAttribute('data-msp-styles', '');
    if (isDocument(root)) root.head.appendChild(style);
    else root.prepend(style);
    update = () => {
      style.textContent = styleText();
    };
  }
  // Registrations often come in a burst (a view's modules): apply them once.
  let queued = false;
  const schedule = () => {
    if (queued) return;
    queued = true;
    queueMicrotask(() => {
      queued = false;
      update();
    });
  };
  roots.set(root, schedule);
  update();
}

/** Stops updating a shadow root's styles (its viewer is gone). */
export function releaseStyles(root: StyleRoot): void {
  roots.delete(root);
}

function isDocument(root: StyleRoot): root is Document {
  return (root as Document).nodeType === 9;
}
