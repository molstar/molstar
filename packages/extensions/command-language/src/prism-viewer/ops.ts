/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure, StructureElement } from '@molstar/model/model/structure';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type { ApiContext, ApiResult, RunOptions } from '@molstar/prism-viewer/api/area';
import { loadedStructures } from '@molstar/prism-viewer/api/scope';
import type { Sel } from '../language/ir.js';
import { atomCount, lociOf } from './atoms.js';
import { objectsOf, type StateHolder } from './state.js';

// What the host's parts share (_spec §5.1): the run's API (calls join `script.run`'s step), the session's state, and
// a selection's atoms per structure. A failed call stops the command with that call's message and hint.

export class HostFailure extends Error {
  constructor(
    message: string,
    readonly hint?: string,
  ) {
    super(message);
  }
}

export type Ok = Extract<ApiResult, { ok: true }>;

/** A selection's atoms on one structure. */
export interface Where {
  readonly object: string;
  readonly ref: string;
  readonly structure: StructureRef;
  readonly data: Structure;
  readonly loci: StructureElement.Loci;
}

export class Ops {
  constructor(
    readonly ctx: ApiContext,
    readonly holder: StateHolder,
  ) {}

  get plugin() {
    return this.ctx.plugin;
  }

  /** Runs an operation as part of `script.run` (its step); a failure stops the command. */
  async run(id: string, values: Readonly<Record<string, unknown>>, options?: RunOptions): Promise<Ok> {
    if (!this.ctx.api.describe(id)) {
      throw new HostFailure(`This viewer cannot ${id.replace('.', ' ')} (it has no '${id}').`);
    }
    const r = await this.ctx.api.run(id, values, options);
    if (!r.ok) throw new HostFailure(r.message, r.hint);
    return r;
  }

  structures(): readonly StructureRef[] {
    return loadedStructures(this.ctx);
  }

  structure(ref: string): StructureRef | undefined {
    return this.structures().find((s) => s.cell.transform.ref === ref);
  }

  objects() {
    return objectsOf(this.ctx, this.holder.value);
  }

  /** The structure ref of an object. */
  ref(object: string): string {
    const lower = object.toLowerCase();
    const o = this.objects().find((x) => x.name.toLowerCase() === lower);
    if (!o)
      throw new HostFailure(
        `There is no object '${object}'.`,
        `Objects: ${
          this.objects()
            .map((x) => x.name)
            .join(', ') || 'none'
        }.`,
      );
    return o.ref;
  }

  /** A selection's atoms on each structure where it selects some. */
  where(sel: Sel): Where[] {
    const out: Where[] = [];
    const byRef = new Map(this.structures().map((s) => [s.cell.transform.ref, s]));
    for (const o of this.objects()) {
      const s = byRef.get(o.ref);
      const data = s?.cell.obj?.data;
      if (!s || !data) continue;
      const loci = lociOf(sel, o.name, data);
      if (atomCount(loci) > 0) out.push({ object: o.name, ref: o.ref, structure: s, data, loci });
    }
    return out;
  }

  /** As `where`, failing when the selection selects nothing anywhere. */
  some(sel: Sel): Where[] {
    const w = this.where(sel);
    if (!w.length) {
      throw new HostFailure(
        `'${sel.text}' selects nothing${this.structures().length ? '' : ': nothing is loaded'}.`,
        this.structures().length ? 'Check the selection.' : 'Load a structure first: fetch 1cbs.',
      );
    }
    return w;
  }
}
