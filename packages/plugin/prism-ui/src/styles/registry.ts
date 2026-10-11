/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
const documents = new Map<Document, () => void>();

/** Registers (or replaces) the CSS of a file; adopted documents update. */
export function defineStyles(id: string, css: string): void {
  if (entries.get(id) === css) return;
  entries.set(id, css);
  for (const update of documents.values()) update();
}

/** All registered CSS, in registration order, after the layer order. */
export function styleText(): string {
  return [LayerOrder, ...entries.values()].join('\n');
}

/**
 * Keeps a document's styles in sync with the registry (idempotent). Uses a constructed stylesheet where supported,
 * otherwise a `<style>` element in the head.
 */
export function adoptStyles(doc: Document = document): void {
  if (documents.has(doc)) return;
  const win = doc.defaultView;
  let update: () => void;
  if (win && 'adoptedStyleSheets' in doc && typeof win.CSSStyleSheet?.prototype.replaceSync === 'function') {
    const sheet = new win.CSSStyleSheet();
    doc.adoptedStyleSheets = [...doc.adoptedStyleSheets, sheet];
    update = () => sheet.replaceSync(styleText());
  } else {
    const style = doc.createElement('style');
    style.setAttribute('data-msp-styles', '');
    doc.head.appendChild(style);
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
  documents.set(doc, schedule);
  update();
}
