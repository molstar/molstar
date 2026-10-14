/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `TreeSchema` with MVS's names (generic/tree-schema.ts there): a root kind and, per node kind, `description`,
// `parent` and `params`. Added: `ref` (`implied` for nodes whose ref is their kind), `count` under the parent,
// `exclusive` param sets, the schema of a param that a registry entry owns, and a node's own rules.

import type { DocNode } from "../doc.js";
import type { SvsRegistries } from "../registry.js";
import type { CheckContext } from "./field-schema.js";
import type { ParamsSchema } from "./params-schema.js";

export interface IndexedNode {
  readonly node: DocNode;
  /** Its ref; "" when it has none (reported). */
  readonly ref: string;
  readonly kind: string;
  readonly params: Readonly<Record<string, unknown>>;
  /** The parent's ref; null for the root. */
  readonly parent: string | null;
  readonly path: string;
}

export interface TreeContext extends CheckContext {
  /** Every node with a ref, by ref. */
  readonly nodes: ReadonlyMap<string, IndexedNode>;
}

export interface NodeSchema {
  readonly description: string;
  /** Kinds this kind may sit under; empty: the root only. */
  readonly parent: readonly string[];
  readonly params: ParamsSchema;
  /** `implied`: the ref may be left out and is then the kind (`root`, `data`); it takes no other. */
  readonly ref?: "required" | "implied";
  /** How many under one parent: exactly one, at most one, or one or more. */
  readonly count?: "one" | "at_most_one" | "some";
  /** Sets of params of which exactly one is given (not null). */
  readonly exclusive?: readonly (readonly string[])[];
  /** A param a registry entry owns (`options`), and its schema for these params; undefined: not checked. */
  readonly owned?: {
    readonly param: string;
    schema(params: Record<string, unknown>, reg: SvsRegistries): ParamsSchema | undefined;
  };
  /** The kind's own rules, run after its params are checked. */
  check?(n: IndexedNode, ctx: TreeContext): void;
}

export type NodeSchemas = { readonly [kind: string]: NodeSchema };
export interface TreeSchema<N extends NodeSchemas = NodeSchemas> {
  readonly rootKind: string;
  readonly nodes: N;
}
export function TreeSchema<N extends NodeSchemas>(schema: { rootKind: keyof N & string; nodes: N }): TreeSchema<N> {
  return schema;
}

/** A node's children as a list: anything but an array is none (reported once, where the node is indexed). */
export const nodeList = (v: unknown): DocNode[] => (Array.isArray(v) ? v : []);
