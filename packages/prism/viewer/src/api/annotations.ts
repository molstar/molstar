/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Structure } from '@molstar/model/model/structure';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { type ApiContext, ApiArea, type ApiResult, ApiResult as Result } from './area.js';
import { isFailure, labelOf, loadedStructures, someStructures } from './scope.js';

// `annotations` (operations.md §2; README decision 27): what a structure carries and the viewer can show of it: quality
// scores, validation, B-factors, secondary structure. Each kind is a provider; extensions register theirs
// (`ViewerExtension.annotations`), so v5's "Automatic (w/ Annotation)" preset becomes a command, not a copy.

export interface AnnotationProvider {
  /** `b-factor`, `plddt`, `validation`, … */
  readonly id: string;
  readonly title: string;
  readonly description?: string;
  /** Whether the structure carries it. */
  isApplicable(structure: Structure): boolean;
  /** Shows it on the structure (through other operations, so it is part of this step). */
  apply(ctx: ApiContext, structure: string): Promise<ApiResult>;
}

/** The viewer's own: what a coordinate file carries. */
export const BuiltInAnnotations: readonly AnnotationProvider[] = [
  {
    id: 'b-factor',
    title: 'B-factors',
    description: 'Atomic displacement (uncertainty), as color',
    isApplicable: (s) => {
      const col = s.model.atomicConformation.B_iso_or_equiv;
      return !!col && col.isDefined && col.rowCount > 0 && col.value(0) !== col.value(col.rowCount - 1);
    },
    apply: (ctx, structure) => ctx.api.representation.color({ target: 'all', theme: 'uncertainty', structure }),
  },
  {
    id: 'secondary-structure',
    title: 'Secondary structure',
    description: 'Helices and strands, as color',
    isApplicable: (s) => s.polymerResidueCount > 0,
    apply: (ctx, structure) =>
      ctx.api.representation.color({ target: 'polymer', theme: 'secondary-structure', structure }),
  },
];

export class AnnotationsApi extends ApiArea {
  readonly id = 'annotations';
  readonly #providers = observableState<readonly AnnotationProvider[]>(BuiltInAnnotations);

  /** The providers registered, the viewer's first. */
  get providers(): ReadonlyState<readonly AnnotationProvider[]> {
    return this.#providers.asReadonly;
  }

  /** Adds a provider (an extension's); returns its removal. An id taken throws. */
  register(provider: AnnotationProvider): () => void {
    if (this.#providers.value.some((p) => p.id === provider.id))
      throw new Error(`An annotation '${provider.id}' exists already`);
    this.#providers.set([...this.#providers.value, provider]);
    return () => this.#providers.set(this.#providers.value.filter((p) => p !== provider));
  }

  readonly list = this.operation({
    id: 'annotations.list',
    title: 'List what the structures carry',
    params: { structure: PD.Text('', { description: 'A structure by label or ref; empty: every one' }) },
    read: true,
    run: async ({ structure }, ctx) => {
      const list = structure ? someStructures(ctx, structure) : loadedStructures(ctx);
      if (isFailure(list)) return list;
      const value = list.map((s) => ({
        structure: labelOf(ctx, s),
        annotations: this.#providers.value
          .filter((p) => s.cell.obj?.data && p.isApplicable(s.cell.obj.data))
          .map((p) => ({ id: p.id, title: p.title })),
      }));
      const text = value
        .map(
          (v) => `${v.structure}: ${v.annotations.map((a) => `${a.title} (${a.id})`).join(', ') || 'nothing to show'}`,
        )
        .join('\n');
      return Result.ok(text || 'Nothing is loaded.', {}, { value });
    },
  });

  readonly apply = this.operation({
    id: 'annotations.apply',
    title: 'Show what a structure carries',
    params: {
      annotation: PD.Text('', { description: 'An annotation by id: b-factor, secondary-structure, …' }),
      structure: PD.Text('', { description: 'A structure by label or ref; empty: every one that carries it' }),
    },
    run: async ({ annotation, structure }, ctx) => {
      const provider = this.#providers.value.find((p) => p.id === annotation);
      if (!provider) {
        return Result.fail(
          'not-found',
          `There is no annotation '${annotation}'.`,
          `Annotations: ${this.#providers.value.map((p) => p.id).join(', ')}.`,
        );
      }
      const list = someStructures(ctx, structure || undefined);
      if (isFailure(list)) return list;
      const carrying = list.filter((s) => s.cell.obj?.data && provider.isApplicable(s.cell.obj.data));
      if (!carrying.length) {
        return Result.fail(
          'not-found',
          `${list.map((s) => labelOf(ctx, s)).join(', ')} carries no ${provider.title.toLowerCase()}.`,
        );
      }
      const nodes: string[] = [];
      for (const s of carrying) {
        const r = await provider.apply(ctx, s.cell.transform.ref);
        if (!r.ok) return r;
        nodes.push(...(r.changed.nodes ?? []));
      }
      return Result.ok(
        `Showing ${provider.title.toLowerCase()} on ${carrying.map((s) => labelOf(ctx, s)).join(', ')}.`,
        { nodes },
      );
    },
  });
}
