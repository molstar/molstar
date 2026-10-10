/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expressionValidationIssues } from '@molstar/query-language/language/validation';
import type { ExpressionValidationOptions } from '@molstar/query-language/language/validation';
import type { Tree } from './tree/generic/tree-schema.js';

/** Validate MolQL in scene/animation trees, including primitive positions, without loading the query runtime. */
export function molQLValidationIssues(
  tree: Tree,
  options: ExpressionValidationOptions = {},
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
    visitParams(node.params, `${path}.params`);
    node.children?.forEach((child, i) => visitNode(child, `${path}.children[${i}]`));
  };
  visitNode(tree, path);
  return issues.length ? issues : undefined;
}
