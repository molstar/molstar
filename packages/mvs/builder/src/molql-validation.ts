/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expressionValidationIssues } from '@molstar/query-language/language/validation';
import type { ExpressionValidationOptions } from '@molstar/query-language/language/validation';
import type { Tree, TreeSchema } from './tree/generic/tree-schema.js';
import * as iots from 'io-ts';
import { PrimitiveFieldTypes } from './tree/generic/field-schema.js';
import { MolQLExpressionT, PrimitiveMolQLExpressionT } from './tree/mvs/param-types.js';

export interface MolQLValidationOptions extends ExpressionValidationOptions {
  /** Only inspect MolQL fields declared by this schema, including nested field types. */
  schema?: TreeSchema;
}

/** Validate MolQL in scene/animation trees, including primitive positions, without loading the query runtime. */
export function molQLValidationIssues(
  tree: Tree,
  options: MolQLValidationOptions = {},
  path = 'tree',
): string[] | undefined {
  const issues: string[] = [];
  visitMolQLExpressions(
    tree,
    (expression, path) => {
      const result = expressionValidationIssues(expression, options);
      if (result) issues.push(...result.map((issue) => `${path}${issue.slice('expression'.length)}`));
    },
    options,
    path,
  );
  return issues.length ? issues : undefined;
}

/** Visit expression values and their paths. With a schema, ignored fields are never traversed.
 * Call after structural validation; without a schema, discover wrappers in arbitrary tree parameters. */
export function visitMolQLExpressions(
  tree: Tree,
  visit: (expression: unknown, path: string) => void,
  options: { schema?: TreeSchema } = {},
  path = 'tree',
): void {
  const visitParams = (value: unknown, path: string, type?: iots.Mixed) => {
    // (only objects can hold an expression: primitives return before any type check, which matters for large color
    // maps, whose values would each be matched against the many cases of a color type, and coordinate lists)
    if (value === null || typeof value !== 'object') return;
    if (type) {
      // (a type that cannot hold an expression is not walked: no `is` checks of large values against it)
      if (!mayContainMolQL(type)) return;
      if (type === MolQLExpressionT || type === PrimitiveMolQLExpressionT) {
        if (value && typeof value === 'object' && Object.hasOwn(value, 'molql')) {
          visit((value as { molql: unknown }).molql, `${path}.molql`);
        }
      } else if (type instanceof iots.UnionType) {
        const branch = type.types.find((branch: iots.Mixed) => branch.is(value));
        if (branch) visitParams(value, path, branch);
      } else if (type instanceof iots.IntersectionType) {
        for (const branch of type.types) visitParams(value, path, branch);
      } else if (Array.isArray(value)) {
        if (type instanceof iots.ArrayType) {
          value.forEach((item, i) => visitParams(item, `${path}[${i}]`, type.type));
        } else if (type instanceof iots.TupleType) {
          type.types.forEach((itemType: iots.Mixed, i: number) => visitParams(value[i], `${path}[${i}]`, itemType));
        }
      } else if (value && typeof value === 'object') {
        const object = value as Record<string, unknown>;
        if (type instanceof iots.InterfaceType || type instanceof iots.PartialType) {
          for (const key of Object.keys(type.props)) visitParams(object[key], `${path}.${key}`, type.props[key]);
        } else if (type instanceof iots.DictionaryType) {
          for (const key of Object.keys(object)) visitParams(object[key], `${path}.${key}`, type.codomain);
        }
      }
      return;
    }
    if (Array.isArray(value)) {
      for (let i = 0; i < value.length; i++) {
        const item = value[i];
        if (item && typeof item === 'object') visitParams(item, `${path}[${i}]`);
      }
    } else if (value && typeof value === 'object') {
      const object = value as Record<string, unknown>;
      if (
        Object.hasOwn(object, 'molql') &&
        Object.keys(object).every((key) => key === 'molql' || key === 'structure_ref')
      ) {
        visit(object.molql, `${path}.molql`);
        return;
      }
      for (const key of Object.keys(object)) {
        const item = object[key];
        if (item && typeof item === 'object') visitParams(item, `${path}.${key}`);
      }
    }
  };
  const visitNode = (node: Tree, path: string) => {
    const params = (node.params ?? {}) as Record<string, unknown>;
    const schema = options.schema?.nodes[node.kind]?.params;
    const fields =
      schema?.type === 'simple' ? schema.fields : schema?.cases[params[schema.discriminator] as string]?.fields;
    const keys = options.schema ? Object.keys(fields ?? {}) : Object.keys(params);
    for (const key of keys) visitParams(params[key], `${path}.params.${key}`, fields?.[key]?.type);
    node.children?.forEach((child, i) => visitNode(child, `${path}.children[${i}]`));
  };
  visitNode(tree, path);
}

const MayContainMolQL = new WeakMap<iots.Mixed, boolean>();

/** Whether a value of `type` can hold a MolQL expression; conservative (true) for types it does not know. */
function mayContainMolQL(type: iots.Mixed): boolean {
  const cached = MayContainMolQL.get(type);
  if (cached !== undefined) return cached;
  MayContainMolQL.set(type, true); // (a recursive type, while it is being decided)
  let result: boolean;
  if (type === MolQLExpressionT || type === PrimitiveMolQLExpressionT) {
    result = true;
  } else if (
    type instanceof iots.UnionType ||
    type instanceof iots.IntersectionType ||
    type instanceof iots.TupleType
  ) {
    result = (type.types as iots.Mixed[]).some(mayContainMolQL);
  } else if (type instanceof iots.ArrayType) {
    result = mayContainMolQL(type.type);
  } else if (type instanceof iots.DictionaryType) {
    result = mayContainMolQL(type.codomain);
  } else if (type instanceof iots.InterfaceType || type instanceof iots.PartialType) {
    result = Object.values(type.props as Record<string, iots.Mixed>).some(mayContainMolQL);
  } else if (
    PrimitiveFieldTypes.has(type) ||
    type instanceof iots.StringType ||
    type instanceof iots.NumberType ||
    type instanceof iots.BooleanType ||
    type instanceof iots.NullType ||
    type instanceof iots.UndefinedType ||
    type instanceof iots.LiteralType ||
    type instanceof iots.KeyofType
  ) {
    result = false;
  } else if (type instanceof iots.RefinementType) {
    result = mayContainMolQL(type.type);
  } else {
    result = true;
  }
  MayContainMolQL.set(type, result);
  return result;
}
