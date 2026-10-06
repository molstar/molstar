import { Expression } from '@molstar/model/script/language/expression';
import { compile } from '@molstar/model/script/runtime/query/base';
import type { Tree } from '@molstar/mvs-builder/tree/generic/tree-schema';
import type { PluginContext } from '@molstar/plugin/context';

/** Preserve runtime MolScript compiler validation after the builder becomes standalone. */
export function molQLValidationIssues(tree: Tree): string[] | undefined {
  const issues: string[] = [];
  const validationCache = new WeakMap<object, true | string>();
  const visitParams = (value: unknown, path: string) => {
    if (Array.isArray(value)) {
      value.forEach((item, i) => visitParams(item, `${path}[${i}]`));
    } else if (value && typeof value === 'object') {
      const object = value as Record<string, unknown>;
      const keys = Object.keys(object);
      if (
        Object.prototype.hasOwnProperty.call(object, 'molql') &&
        keys.every((key) => key === 'molql' || key === 'structure_ref')
      ) {
        const expression = object.molql;
        if (!Expression.is(expression) || !Expression.isApply(expression)) {
          issues.push(`${path}.molql: MolQL expression must be a MolScript application`);
        } else {
          const cacheKey = expression !== null && typeof expression === 'object' ? expression : undefined;
          let result = cacheKey ? validationCache.get(cacheKey) : undefined;
          if (result === undefined) {
            try {
              compile(expression);
              result = true;
            } catch (e) {
              result = `Invalid MolQL expression: ${e instanceof Error ? e.message : String(e)}`;
            }
            if (cacheKey) validationCache.set(cacheKey, result);
          }
          if (result !== true) issues.push(`${path}.molql: ${result}`);
        }
      }
      for (const key of Object.keys(object)) visitParams(object[key], `${path}.${key}`);
    }
  };
  const visitNode = (node: Tree, path: string) => {
    visitParams(node.params, `${path}.params`);
    node.children?.forEach((child, i) => visitNode(child, `${path}.children[${i}]`));
  };
  visitNode(tree, 'tree');
  return issues.length ? issues : undefined;
}

export function validateMolQLTree(tree: Tree, label: string, plugin: PluginContext): void {
  const issues = molQLValidationIssues(tree);
  if (!issues) return;
  console.warn(`Invalid ${label} tree:\n${issues.join('\n')}`);
  console.error(`${label} tree validation issues:`);
  plugin.log.error(`${label} tree validation issues:`);
  for (const line of issues) {
    console.error(' ', line);
    plugin.log.error(line);
  }
  throw new Error('FormatError');
}
