/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expressionValidationIssues } from '@molstar/query-language/language/validation';
import type { ExpressionValidationOptions } from '@molstar/query-language/language/validation';
import type { Tree, TreeSchema } from './tree/generic/tree-schema.js';

export interface MolQLValidationOptions extends ExpressionValidationOptions {
  /** Only inspect parameters declared by this schema. Extra parameters are handled by structural validation. */
  schema?: TreeSchema;
}

/** Validate MolQL in scene/animation trees, including primitive positions, without loading the query runtime. */
export function molQLValidationIssues(
  tree: Tree,
  options: MolQLValidationOptions = {},
  path = 'tree',
): string[] | undefined {
  const issues: string[] = [];
  const visitParams = (value: unknown, path: string) => {
    if (Array.isArray(value)) {
      value.forEach((item, i) => visitParams(item, `${path}[${i}]`));
    } else if (value && typeof value === 'object') {
      const object = value as Record<string, unknown>;
      if (
        Object.hasOwn(object, 'molql') &&
        Object.keys(object).every((key) => key === 'molql' || key === 'structure_ref')
      ) {
        const result = expressionValidationIssues(object.molql, options);
        if (result) issues.push(...result.map((issue) => `${path}.molql${issue.slice('expression'.length)}`));
        return;
      }
      for (const key of Object.keys(object)) visitParams(object[key], `${path}.${key}`);
    }
  };
  const visitNode = (node: Tree, path: string) => {
    const params = (node.params ?? {}) as Record<string, unknown>;
    const schema = options.schema?.nodes[node.kind]?.params;
    const fields =
      schema?.type === 'simple' ? schema.fields : schema?.cases[params[schema.discriminator] as string]?.fields;
    const keys = options.schema ? Object.keys(fields ?? {}) : Object.keys(params);
    for (const key of keys) visitParams(params[key], `${path}.params.${key}`);
    node.children?.forEach((child, i) => visitNode(child, `${path}.children[${i}]`));
  };
  visitNode(tree, path);
  return issues.length ? issues : undefined;
}
