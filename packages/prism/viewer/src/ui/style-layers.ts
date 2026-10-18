/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { mapState } from '@molstar/prism-ui/reactive/map-state';
import { type ObservableState, observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import type { ApiResult, ViewerApi } from '../api/area.js';
import type { StyleLayerItem, StyleLayersModel, StyleLayerView } from '../model/style-layers.js';
import { sameItems } from '../sources/state.js';
import { type StyleLayerInfo, type StyleLayerKind, styleLayerId, styleLayerKey } from '../sources/style-layers.js';
import { Debouncer } from './debounce.js';

// The model of the applied style layers under the Selection tab's style form (state.md: views change the plugin only
// through operations): the rows follow `sources.layers` (watched: nothing listens while the section is folded); an
// edit runs `representation.layer-update` (debounced, one step per settled edit), remove `layer-remove`, Clear all
// `representation.reset` with every layer. DOM-free.

const Icons: Record<StyleLayerKind, IconName> = {
  color: 'paintbrush',
  transparency: 'eye-off',
  emissive: 'sun',
  substance: 'sparkles',
  clipping: 'crop',
};

const Titles: Record<StyleLayerKind, string> = {
  color: 'Color',
  transparency: 'Transparency',
  emissive: 'Glow',
  substance: 'Material',
  clipping: 'Excluded from clipping',
};

/** The row's value, short. */
function valueOf(l: StyleLayerInfo): string {
  switch (l.kind) {
    case 'color':
      return l.clear ? 'color reset' : '';
    case 'transparency':
      return `${Math.round((l.amount ?? 0) * 100)}% transparent`;
    case 'emissive':
      return `glow ${l.amount ?? 0}`;
    case 'substance':
      return l.clear ? 'material reset' : 'material';
    case 'clipping':
      return 'unclipped';
  }
}

function viewOf(l: StyleLayerInfo, several: boolean): StyleLayerView {
  const head = l.clear ? `${Titles[l.kind]} reset` : `${Titles[l.kind]} ${l.value}`;
  return {
    ...(l.color ? { swatch: l.color } : {}),
    icon: l.clear ? 'rotate-ccw' : Icons[l.kind],
    where: several ? `${l.structure.label} · ${l.where}` : l.where,
    value: valueOf(l),
    detail: `${head} on ${l.line} of ${l.structure.label}`,
  };
}

const sameView = (a: StyleLayerView, b: StyleLayerView) =>
  a.swatch === b.swatch && a.icon === b.icon && a.where === b.where && a.value === b.value && a.detail === b.detail;

/** The value an edit sends (`representation.layer-update`'s `value`), and the key it gives the layer. */
function sent(kind: StyleLayerKind, values: { value: unknown }): { value: string; key: string } {
  if (kind === 'color') {
    const color = values.value as Color;
    return { value: Color.toHexStyle(color), key: styleLayerKey(kind, { color }) };
  }
  const n = Number(values.value);
  return { value: String(n), key: styleLayerKey(kind, { value: n }) };
}

class Row implements StyleLayerItem {
  readonly #view: ObservableState<StyleLayerView>;
  #info: StyleLayerInfo;
  #form: ParamsFormModel<any> | undefined;
  /** The edits in flight, one after the other (each names the layer by the id the one before gave it). */
  chain: Promise<void> = Promise.resolve();

  constructor(
    readonly key: string,
    info: StyleLayerInfo,
    several: boolean,
    private readonly owner: { edited(row: Row): void; remove(row: Row): void },
  ) {
    this.#info = info;
    this.#view = observableState(viewOf(info, several), { compare: sameView });
  }

  /** The layer's id now (an edit moves it on). */
  id: string = '';

  get info(): StyleLayerInfo {
    return this.#info;
  }

  get view(): ReadonlyState<StyleLayerView> {
    return this.#view.asReadonly;
  }

  set(info: StyleLayerInfo, several: boolean): void {
    this.#info = info;
    this.id = info.id;
    this.#view.set(viewOf(info, several));
  }

  get editable(): boolean {
    const l = this.#info;
    return !l.clear && (l.kind === 'color' || l.kind === 'transparency' || l.kind === 'emissive');
  }

  editor(): ParamsFormModel<any> | undefined {
    const l = this.#info;
    if (!this.editable) return undefined;
    const current = l.kind === 'color' ? Color.fromHexStyle(l.color ?? '#000000') : (l.amount ?? 0);
    if (!this.#form) {
      const param =
        l.kind === 'color'
          ? PD.Color(current as Color, { label: 'Color' })
          : PD.Numeric(
              current as number,
              { min: 0, max: 1, step: 0.05 },
              { label: l.kind === 'transparency' ? 'Transparency' : 'Glow' },
            );
      this.#form = new ParamsFormModel({ value: param }, undefined, { onChange: () => this.owner.edited(this) });
    } else {
      // Opened again: the value the state has now (an Undo may have changed it).
      this.#form.load({ value: current });
    }
    return this.#form;
  }

  /** The value the open editor holds. */
  get edit(): { value: unknown } | undefined {
    return this.#form?.values.value as { value: unknown } | undefined;
  }

  remove(): void {
    this.owner.remove(this);
  }
}

export interface StyleLayersHost {
  readonly model: StyleLayersModel;
  dispose(): void;
}

/**
 * The applied layers over `layers` (the viewer's source), changed through `api` (bound to the person). `report`
 * says a failed run where the person looks.
 */
export function createStyleLayersModel(
  api: ViewerApi,
  layers: ReadonlyState<readonly StyleLayerInfo[]>,
  report: (r: ApiResult) => void = () => {},
): StyleLayersHost {
  // The rows by the layer's id now; an edit in flight files its row under the id it will give (`pending`), so the read
  // that brings the new id keeps the row (and its open editor) instead of making a new one.
  const rows = new Map<string, Row>();
  const pending = new Map<string, Row>();

  const send = async (row: Row) => {
    const edit = row.edit;
    if (!edit) return;
    const from = row.id;
    const kind = row.info.kind;
    const { value, key } = sent(kind, edit);
    const to = styleLayerId(kind, row.info.structure.ref, key);
    if (to === from) return;
    pending.set(to, row);
    const r = await api.representation.layerUpdate({ layer: from, value });
    if (!r.ok) {
      pending.delete(to);
      report(r);
      return;
    }
    // Landed: the row is the new id's, whether or not the list has been read since.
    if (rows.get(from) === row) rows.delete(from);
    rows.set(to, row);
    row.id = to;
  };
  const debouncer = new Debouncer<Row>((row) => {
    row.chain = row.chain.then(() => send(row));
  });

  const owner = {
    edited: (row: Row) => debouncer.schedule(row),
    remove: (row: Row) => {
      debouncer.cancel(row);
      void row.chain.then(async () => {
        const r = await api.representation.layerRemove({ layer: row.id });
        if (!r.ok) report(r);
      });
    },
  };

  const rowsOf = (list: readonly StyleLayerInfo[]): readonly StyleLayerItem[] => {
    const several = new Set(list.map((l) => l.structure.ref)).size > 1;
    const live = new Set<string>();
    const out: Row[] = [];
    for (const info of list) {
      let row = rows.get(info.id);
      if (!row) {
        row = pending.get(info.id);
        if (row) {
          pending.delete(info.id);
          if (rows.get(row.id) === row) rows.delete(row.id);
        } else {
          row = new Row(info.id, info, several, owner);
        }
        rows.set(info.id, row);
      }
      row.set(info, several);
      live.add(info.id);
      out.push(row);
    }
    if (rows.size > live.size) {
      const waiting = new Set(pending.values());
      for (const [id, row] of rows) if (!live.has(id) && !waiting.has(row)) rows.delete(id);
    }
    return out;
  };

  const items = mapState(layers, rowsOf, (a, b) => sameItems(a, b));
  return {
    model: {
      items,
      clear: () => {
        debouncer.dispose();
        void api.representation.reset({ layer: 'all' }).then((r) => {
          if (!r.ok) report(r);
        });
      },
    },
    dispose: () => debouncer.dispose(),
  };
}
