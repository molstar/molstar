/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Callable symbol definitions (seqspec:seqql#overview, seqspec:seqql#extension): a symbol called with its arguments
// builds the canonical JSON application, as MolQL's `MSymbol` does; defaults are omitted, so equal queries give equal
// JSON (seqspec:seqql#names).
import type { Apply, Expression, Keywords, SymbolDef, SymbolSpec } from "./types.js";

type ArgValue = Expression | readonly Expression[];

function isKeywords(x: unknown): x is Keywords {
  return typeof x === "object" && x !== null && !Array.isArray(x) && !("head" in x);
}

function apply(def: SymbolSpec, input: readonly (Expression | Keywords)[]): Apply {
  const kw = input.length > 0 && isKeywords(input[input.length - 1]) ? (input[input.length - 1] as Keywords) : {};
  const pos = input.length > 0 && isKeywords(input[input.length - 1]) ? input.slice(0, -1) : input;
  const args: Record<string, ArgValue> = {};
  let next = 0;
  let index = 0;
  for (const a of def.args) {
    if (!a.positional) continue;
    const key = String(index++);
    if (a.rest) {
      if (next < pos.length) args[key] = pos.slice(next) as Expression[];
      next = pos.length;
    } else if (next < pos.length) {
      const v = pos[next++] as Expression;
      if (v !== a.default) args[key] = v;
    }
  }
  if (next < pos.length || pos.some(isKeywords)) {
    throw new Error(`${def.name}: too many positional arguments; ${def.signature}`);
  }
  for (const key of Object.keys(kw)) {
    if (!def.args.some((a) => !a.positional && a.name === key)) throw new Error(`${def.name}: no keyword :${key}`);
  }
  for (const a of def.args) {
    const v = kw[a.name];
    if (!a.positional && v !== undefined && v !== a.default) args[a.name] = v;
  }
  return Object.keys(args).length > 0 ? { head: def.name, args } : { head: def.name };
}

/** A symbol from its definition: the spec's fields on a function that builds its applications. */
export function defineSymbol<const N extends string>(
  spec: SymbolSpec & { readonly name: N },
): SymbolDef & { readonly name: N } {
  const call = (...args: readonly (Expression | Keywords)[]): Apply =>
    spec.normalize ? spec.normalize(apply(spec, args)) : apply(spec, args);
  // A function's own `name` is read-only: redefined, not assigned.
  const { name, ...rest } = spec;
  Object.defineProperty(call, "name", { value: name, enumerable: true });
  return Object.freeze(Object.assign(call, rest)) as never;
}
