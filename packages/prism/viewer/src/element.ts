/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { UITheme } from '@molstar/prism-ui/shell/theme';
import { DefaultPrismViewerSpec } from './default-spec.js';
import { type CreatePrismViewerOptions, createPrismViewer } from './page.js';
import type { PrismViewerSpec } from './spec.js';
import type { PrismViewer } from './viewer.js';

// `<prism-viewer>` (app.md §3): the element wraps `createPrismViewer`, the element itself its target (it gets the
// shadow root the shell mounts in, so host CSS stops there and the `--msp-*` tokens cross it). Its `spec` property takes
// a `PrismViewerSpec`; the attributes cover the common options: `pdb` (an entry to load; a new value adds that one),
// `theme` and `hide-controls` (read when the viewer is made). Several on a page each own their viewer; one leaving the
// page disposes its viewer, and coming back makes a new one. The class is made by a function, so importing this module
// needs no DOM (the headless viewer and Node import the package).

/** The attributes the element reads. */
export const PrismViewerElementAttributes = ['pdb', 'theme', 'hide-controls'] as const;

const Themes: readonly UITheme[] = ['light', 'dark', 'blue', 'contrast', 'auto'];

export interface PrismViewerElementOptions {
  /** The spec while the `spec` property is not set (default: `DefaultPrismViewerSpec()`, the built-ins only). */
  readonly defaultSpec?: () => PrismViewerSpec;
  /** Makes the viewer (default: `createPrismViewer`); a host or a test passes its own. */
  readonly create?: (o: CreatePrismViewerOptions) => Promise<PrismViewer>;
}

export interface PrismViewerElement extends HTMLElement {
  /** The spec the viewer is made from; setting it while connected makes the viewer again. */
  spec: PrismViewerSpec | undefined;
  /** The viewer, while connected and made. */
  readonly viewer: PrismViewer | undefined;
  /**
   * The viewer once it is made (the current one; a new promise after the element comes back or the spec changes);
   * rejects when making it fails (as the `error` event).
   */
  readonly ready: Promise<PrismViewer>;
}

/**
 * The spec with the attributes' options over it: `hide-controls` starts the rail collapsed (as the URL's), `theme` is
 * the UI theme. The given spec's own values stay where the attributes say nothing.
 */
export function specWithAttributes(
  spec: PrismViewerSpec,
  attributes: { readonly theme?: string | null; readonly hideControls?: string | null },
): PrismViewerSpec {
  const theme = Themes.find((t) => t === attributes.theme);
  const hide = attributes.hideControls;
  const collapsed = hide === null || hide === undefined ? undefined : hide !== 'false' && hide !== '0';
  return {
    ...spec,
    ...(collapsed !== undefined ? { options: { ...spec.options, ui: { ...spec.options?.ui, collapsed } } } : {}),
    ...(theme ? { ui: { ...spec.ui, theme } } : {}),
  };
}

/** The element's class, over `options` (made when called, so the module imports without a DOM). */
export function createPrismViewerElementClass(options: PrismViewerElementOptions = {}): new () => PrismViewerElement {
  const create = options.create ?? createPrismViewer;
  const defaultSpec = options.defaultSpec ?? DefaultPrismViewerSpec;

  return class extends HTMLElement implements PrismViewerElement {
    static get observedAttributes(): readonly string[] {
      return PrismViewerElementAttributes;
    }

    #spec: PrismViewerSpec | undefined;
    #viewer: PrismViewer | undefined;
    /** Bumped when the spec changes and on leaving: a viewer made for an older one is disposed as it arrives. */
    #generation = 0;
    /** The generation the current viewer was made for. */
    #made = 0;
    /**
     * Makes and disposals, one at a time and a microtask after the change that asked for them: a spec set in the same
     * task as connecting is the one made, and a viewer is disposed (its persist key freed) before the next is made.
     */
    #queue: Promise<void> = Promise.resolve();
    #ready!: Promise<PrismViewer>;
    #resolve!: (v: PrismViewer) => void;
    #reject!: (e: unknown) => void;
    #settled = false;
    #loaded: string | undefined;

    constructor() {
      super();
      this.#newReady();
    }

    get spec(): PrismViewerSpec | undefined {
      return this.#spec;
    }
    set spec(spec: PrismViewerSpec | undefined) {
      this.#spec = spec;
      this.#generation++;
      this.#renewReady();
      if (this.isConnected) this.#sync();
    }

    get viewer(): PrismViewer | undefined {
      return this.#viewer;
    }

    get ready(): Promise<PrismViewer> {
      return this.#ready;
    }

    connectedCallback(): void {
      this.#sync();
    }

    disconnectedCallback(): void {
      // Moved within the page (disconnected and connected in one task): keep the viewer.
      queueMicrotask(() => {
        if (this.isConnected) return;
        this.#generation++;
        this.#sync();
      });
    }

    attributeChangedCallback(name: string, old: string | null, value: string | null): void {
      if (name === 'pdb' && old !== value && this.#viewer) void this.#load(this.#viewer);
    }

    /** A pending `ready` for the viewer to come (one read before stays with the viewer it was for). */
    #renewReady() {
      if (this.#settled) this.#newReady();
    }

    #newReady() {
      this.#settled = false;
      this.#ready = new Promise((resolve, reject) => {
        this.#resolve = resolve;
        this.#reject = reject;
      });
      // Handled here too: a page that never reads `ready` gets the `error` event, not an unhandled rejection.
      this.#ready.catch(() => {});
    }

    /** Brings the viewer in line with the element: none while it is off the page, one made for the current spec. */
    #sync(): void {
      this.#queue = this.#queue
        .then(async () => {
          const want = this.isConnected;
          if (this.#viewer && (!want || this.#made !== this.#generation)) await this.#dispose();
          if (want && !this.#viewer) await this.#make();
        })
        .catch((e) => console.error(e));
    }

    async #make(): Promise<void> {
      const generation = this.#generation;
      this.#renewReady();
      const shadow = this.shadowRoot ?? this.attachShadow({ mode: 'open' });
      // The shell fills its host (`position: absolute; inset: 0`): the element is a block that holds it.
      if (!shadow.querySelector('style[data-prism-viewer]')) {
        const style = document.createElement('style');
        style.dataset.prismViewer = '';
        style.textContent = ':host { display: block; position: relative; min-height: 320px; overflow: hidden }';
        shadow.prepend(style);
      }
      const spec = specWithAttributes(this.#spec ?? defaultSpec(), {
        theme: this.getAttribute('theme'),
        hideControls: this.getAttribute('hide-controls'),
      });
      let viewer: PrismViewer;
      try {
        viewer = await create({ target: this, spec });
      } catch (error) {
        if (generation === this.#generation) {
          this.#settled = true;
          this.#reject(error);
          this.dispatchEvent(new CustomEvent('error', { detail: { error } }));
          console.error(error);
        }
        return;
      }
      if (generation !== this.#generation) {
        await viewer.dispose();
        return;
      }
      this.#viewer = viewer;
      this.#made = generation;
      this.#loaded = undefined;
      this.#settled = true;
      this.#resolve(viewer);
      this.dispatchEvent(new CustomEvent('ready', { detail: { viewer } }));
      void this.#load(viewer);
    }

    async #dispose(): Promise<void> {
      const viewer = this.#viewer;
      this.#viewer = undefined;
      if (!viewer) return;
      this.#renewReady();
      await viewer.dispose();
    }

    /** The `pdb` attribute's entry, through the loaders feature when the spec has it. */
    async #load(viewer: PrismViewer): Promise<void> {
      const id = this.getAttribute('pdb')?.trim();
      if (!id || id === this.#loaded || !viewer.api.describe('load.entry')) return;
      this.#loaded = id;
      const r = await viewer.api.run('load.entry', { source: 'pdb', id });
      if (!r.ok) console.warn(`<prism-viewer pdb="${id}">: ${r.message}`);
    }
  };
}

/**
 * Defines `<prism-viewer>` (or `tag`) once: a second call, or a page that defined it already, keeps the first.
 * Returns whether this call defined it.
 */
export function definePrismViewerElement(options: PrismViewerElementOptions & { readonly tag?: string } = {}): boolean {
  const tag = options.tag ?? 'prism-viewer';
  if (typeof customElements === 'undefined' || customElements.get(tag)) return false;
  customElements.define(tag, createPrismViewerElementClass(options));
  return true;
}
