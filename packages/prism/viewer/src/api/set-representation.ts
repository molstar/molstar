/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ApiContext, type ApiFailure, ApiResult } from './area.js';
import { mergeParams } from './params.js';
import { checkValues } from './run.js';
import { isFailure } from './scope.js';

// "Set a representation type with patched params" (operations.md §2), the one core of `representation.set`,
// `volume.representation` and `particles.representation`: the type exists, applies to the data (the hint lists the types
// that do), the patch is checked against the type's params (`checkValues`), then the representations of that type that
// exist are updated (equality by the type's params, the merge of `mergeParams`), else one is created.

/** What of a representation type's registry the core uses (the graphics `RepresentationRegistry`'s). */
export interface TypeRegistry<D> {
  readonly types: readonly [string, string][];
  has(type: string): boolean;
  get(type: string): {
    readonly label: string;
    isApplicable(data: D): boolean;
    getParams(themes: any, data: D): PD.Params;
  };
  getApplicableTypes(data: D): [string, string][];
}

/** A representation of the scene hierarchy. */
interface RepresentationLike {
  readonly cell: { readonly transform: { readonly ref: string; readonly params?: unknown } };
}

type Values = { type: { name: string; params: Record<string, any> } } & Record<string, any>;

export interface SetRepresentation<D, R extends RepresentationLike> {
  /** Its subject for the person: the component, volume or list. */
  readonly name: string;
  readonly type: string;
  readonly patch: Record<string, unknown> | undefined;
  readonly registry: TypeRegistry<D>;
  /** The theme registries `getParams` takes. */
  readonly themes: unknown;
  readonly data: D | undefined;
  /** The subject's representations: those of the type are updated. */
  readonly existing: readonly R[];
  /** Makes one when there is none of the type; the nodes it made, or why not. */
  create(provider: ReturnType<TypeRegistry<D>['get']>, typeParams: PD.Params): Promise<readonly string[] | ApiFailure>;
}

const typeNameOf = (r: RepresentationLike) => (r.cell.transform.params as Values | undefined)?.type?.name ?? '';

export async function setRepresentation<D, R extends RepresentationLike>(ctx: ApiContext, o: SetRepresentation<D, R>) {
  const { name, type, registry } = o;
  if (!registry.has(type)) {
    return ApiResult.fail(
      'invalid',
      `There is no representation type '${type}'.`,
      `Types: ${registry.types.map((t) => t[0]).join(', ')}.`,
    );
  }
  const data = o.data;
  if (!data) return ApiResult.fail('blocked', `${name} has no data to show.`);
  const provider = registry.get(type);
  if (!provider.isApplicable(data)) {
    return ApiResult.fail(
      'blocked',
      `${provider.label} does not apply to ${name}.`,
      `Applicable: ${registry
        .getApplicableTypes(data)
        .map((t) => t[0])
        .join(', ')}.`,
    );
  }
  const typeParams = provider.getParams(o.themes, data);
  const patch = o.patch ?? {};
  // Unknown names and wrong values: the same check as a call's own values (it names what the type takes).
  const checked = checkValues(typeParams, patch);
  if (isFailure(checked)) return ApiResult.fail(checked.code, `${provider.label}: ${checked.message}`, checked.hint);

  const existing = o.existing.filter((r) => typeNameOf(r) === type);
  if (existing.length) {
    const changing = existing.filter((r) => {
      const old = (r.cell.transform.params as Values).type.params;
      return !PD.areEqual(typeParams, mergeParams(old, patch), old);
    });
    if (!changing.length) return ApiResult.already(`${name} already shows ${provider.label} that way.`);
    await ctx.plugin.dataTransaction(async () => {
      for (const r of changing) {
        const old = r.cell.transform.params as Values;
        await ctx.plugin.state.updateTransform(ctx.plugin.state.data, r.cell.transform.ref, {
          ...old,
          type: { name: old.type.name, params: mergeParams(old.type.params, patch) },
        });
      }
    });
    return ApiResult.ok(`Updated ${provider.label} on ${name}.`, { nodes: changing.map((r) => r.cell.transform.ref) });
  }

  const made = await o.create(provider, typeParams);
  if (isFailure(made)) return made;
  return ApiResult.ok(`Showed ${name} as ${provider.label}.`, { nodes: [...made] });
}

/** Commits an update that applied a representation; the new node, or a failure (the run's checkpoint restores the scene). */
export async function commitRepresentation(
  ctx: ApiContext,
  update: { commit(): Promise<unknown>; selector: { ref: string } },
  failure: string,
): Promise<readonly string[] | ApiFailure> {
  await update.commit();
  const made = ctx.plugin.state.data.cells.get(update.selector.ref);
  // The state keeps the failed cell.
  if (!made || made.status === 'error')
    return ApiResult.fail('failed', `${failure}: ${made?.errorText ?? 'it failed'}.`);
  return [update.selector.ref];
}
