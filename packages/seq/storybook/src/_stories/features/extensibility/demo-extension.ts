/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  defineExtension,
  type Hit,
  type ItemDescription,
  type LayerType,
  type SchemeEntry,
} from '@molstar/seqx-core/extension';
import type { BehaviorProvider, SeqViewerExtension } from '@molstar/seqx-viewer';

// An example extension, written as an outside package writes one (the core's contract from its one entry,
// `@molstar/seqx-core/extension`, and the viewer's descriptor type): a layer type with a plain param, a key it reads on
// a built-in section, a scheme, a command (from its behavior) and a binding. Its id, `demo`, is its namespace.

/** What the lollipops read of their data: each item's extent and key (inline ranges). */
interface LollipopsPrep {
  readonly owner: string;
  readonly items: readonly { readonly key: string | number; readonly start: number; readonly end: number }[];
  readonly head: number;
}

interface LollipopsLay {
  readonly prep: LollipopsPrep;
  /** Each item's head centre (x), its stem's top, and its extent's x in view; NaN where out of view. */
  xs: number[];
  x0: number[];
  x1: number[];
  stemTop: number;
  bottom: number;
  r: number;
  stem: number;
  head: number;
  key: string;
  section: string;
  system: string;
}

/**
 * `demo:lollipops`: a stem and a head at each range's middle (inline ranges, on the section's own system). Its plain
 * param `head` is the head's radius in px (4 when absent); the section's `demo:stems` (`short`, else tall) how far the
 * stems reach. A reader without the extension draws its `fallback` (`ranges`).
 */
export const LollipopsType: LayerType<LollipopsPrep, LollipopsLay> = {
  name: 'demo:lollipops',
  shape: 'ranges',
  fallback: 'ranges',
  prepare(layer, model) {
    const owner = layer.data?.address ?? '';
    const head = layer.params.head;
    return {
      owner,
      items: model.annotations.get(owner)?.items ?? [],
      head: typeof head === 'number' && head > 0 ? head : 4,
    };
  },
  layout(prep, ctx, _data, _overlay, prev) {
    const lay: LollipopsLay =
      prev?.prep === prep
        ? prev
        : {
            prep,
            xs: [],
            x0: [],
            x1: [],
            stemTop: 0,
            bottom: 0,
            r: 0,
            stem: 0,
            head: 0,
            key: '',
            section: '',
            system: '',
          };
    const { box, axis, section, tokens } = ctx;
    lay.xs.length = lay.x0.length = lay.x1.length = 0;
    for (const it of prep.items) {
      const a = axis.x(it.start);
      const b = axis.x(it.end);
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      const visible = Number.isFinite(lo) && hi >= box.x && lo <= box.x + box.width;
      lay.xs.push(visible ? (lo + hi) / 2 : Number.NaN);
      lay.x0.push(Math.max(lo, box.x));
      lay.x1.push(Math.min(hi, box.x + box.width));
    }
    lay.r = Math.min(prep.head, box.height / 3);
    const short = section.params['demo:stems'] === 'short';
    lay.bottom = box.y + box.height;
    lay.stemTop = short ? box.y + box.height / 2 : box.y + lay.r;
    lay.stem = tokens.colors['seq-range'];
    lay.head = tokens.colors['seq-cat-1'];
    lay.key = ctx.key;
    lay.section = section.address;
    lay.system = section.axis.system;
    return lay;
  },
  draw(lay, pass, out) {
    if (pass !== 'data') return;
    const r = lay.r;
    for (const x of lay.xs) {
      if (Number.isNaN(x)) continue;
      out.fill(x - 0.75, lay.stemTop, 1.5, lay.bottom - lay.stemTop, lay.stem);
      out.rect(x - r, lay.stemTop - r, 2 * r, 2 * r, lay.head, 0, 0, r);
    }
  },
  hit(lay, at): Hit | undefined {
    const r = lay.r + 2;
    for (let i = 0; i < lay.xs.length; i++) {
      const x = lay.xs[i];
      if (Number.isNaN(x) || Math.abs(x - at[0]) > r || at[1] < lay.stemTop - r || at[1] > lay.bottom) continue;
      const it = lay.prep.items[i];
      return {
        layer: lay.key,
        section: lay.section,
        system: lay.system,
        element: { owner: lay.prep.owner, key: it.key },
        hull: [it.start, it.end],
      };
    }
    return undefined;
  },
  items(lay): readonly ItemDescription[] {
    return lay.prep.items.flatMap((it, i) =>
      Number.isNaN(lay.xs[i])
        ? []
        : [{ key: it.key, label: String(it.key), start: it.start, end: it.end, x0: lay.x0[i], x1: lay.x1[i] }],
    );
  },
};

/** `demo:warm`: a category a warm color (a categories scheme a built-in layer's theme names). */
export const WarmScheme: SchemeEntry = {
  name: 'demo:warm',
  label: 'Warm',
  subject: 'categories',
  cvdSafe: false,
  dark: 'derive',
  expand: (p) => ({
    kind: 'field',
    params: {
      name: typeof p.field === 'string' ? p.field : 'category',
      palette: { kind: 'categorical', params: { colors: ['#d9480f', '#f59f00', '#c2255c', '#e8590c'] } },
    },
  }),
};

/** The command it adds: shows the sites (focuses their data), from its behavior (which holds the viewer's models). */
export const FocusSites = 'demo:focus-sites';

const DemoBehavior: BehaviorProvider = {
  name: 'demo:behavior',
  create: (ctx) => ({
    commands: [
      {
        id: FocusSites,
        title: 'Show the sites',
        category: 'Demo',
        run: () => {
          ctx.addresses.focus('/data/sites');
          return { message: 'The sites are shown' };
        },
      },
    ],
  }),
};

/** The example extension: everything it ships named `demo:…`, its params checked where it is composed. */
export const DemoExtension: SeqViewerExtension = defineExtension({
  id: 'demo',
  version: '0.1.0',
  requires: '0.5',
  doc: 'Lollipops, a warm scheme, and a command that shows the sites (Alt+double-click).',
  registry: [
    { kind: 'layer-type', type: LollipopsType },
    { kind: 'scheme', scheme: WarmScheme },
  ],
  params: {
    'demo:lollipops': { head: { type: 'number', min: 1, max: 12, doc: 'The heads’ radius, px (4).' } },
    section: { 'demo:stems': { type: 'string', enum: ['tall', 'short'], doc: 'How far the stems reach (tall).' } },
  },
  behaviors: [{ provider: DemoBehavior }],
  bindings: [{ gesture: 'double-click', modifiers: ['alt'], target: 'any', action: FocusSites }],
});
