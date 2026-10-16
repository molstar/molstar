/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

// What a change of settings says in words, for the History and the result line ("Shadow on, Color white"): the card
// of the Canvas and the Settings tab name what they changed the same way (DOM-free, no plugin).

const words = (key: string) => key.replace(/([a-z0-9])([A-Z])/g, '$1 $2').replace(/^./, (c) => c.toUpperCase());

const labelOf = (key: string, p: PD.Any) => p.label ?? words(key);

/** Most changes a line names; the rest are counted. */
const Limit = 4;

const colorWord = (c: number) => {
  const hex = Color.toHexStyle(c as Color);
  return hex === '#ffffff' ? 'white' : hex === '#000000' ? 'black' : hex;
};

const round = (n: number) => String(Math.round(n * 1000) / 1000);

const plain = (x: unknown): x is Record<string, any> =>
  !!x && typeof x === 'object' && !Array.isArray(x) && Object.getPrototypeOf(x) === Object.prototype;

function one(
  key: string,
  p: PD.Any,
  before: any,
  after: any,
  names: Readonly<Record<string, string>>,
  path: string,
): string[] {
  if (PD.isParamEqual(p, before, after)) return [];
  const label = names[path] ?? labelOf(key, p);
  switch (p.type) {
    case 'boolean':
      return [`${label} ${after ? 'on' : 'off'}`];
    case 'number':
      return [`${label} ${round(after)}`];
    case 'color':
      return [`${label} ${colorWord(after)}`];
    case 'select': {
      const option = p.options.find((o) => o[0] === after);
      return [`${label} ${String(option?.[1] ?? after).toLowerCase()}`];
    }
    case 'text':
      return [`${label} "${String(after).slice(0, 24)}"`];
    case 'group':
      return inner(p.params, before, after, names, path);
    case 'mapped': {
      const name = after?.name;
      if (before?.name !== name) {
        const option = p.select.options.find((o) => o[0] === name);
        return [`${label} ${String(option?.[1] ?? name).toLowerCase()}`];
      }
      const params = p.map(name);
      const changes =
        params.type === 'group' ? inner(params.params, before?.params, after?.params, names, `${path}.params`) : [];
      return changes.length ? [`${label}: ${changes.join(', ')}`] : [`${label} changed`];
    }
    default:
      return [`${label} changed`];
  }
}

function inner(
  params: PD.Params,
  before: unknown,
  after: unknown,
  names: Readonly<Record<string, string>>,
  path: string,
): string[] {
  const a = plain(before) ? before : {};
  const b = plain(after) ? after : {};
  return Object.keys(params).flatMap((key) =>
    one(key, params[key], a[key], b[key], names, path ? `${path}.${key}` : key),
  );
}

/**
 * The fields that differ between two values of `params`, in words ("Shadow on", "Color white", "Camera mode
 * orthographic"); `names` renames a field by its path (`renderer.backgroundColor`). Long lists end in a count.
 */
export function changeWords(
  params: PD.Params,
  before: Record<string, unknown>,
  after: Record<string, unknown>,
  names: Readonly<Record<string, string>> = {},
): string[] {
  const all = inner(params, before, after, names, '');
  return all.length > Limit ? [...all.slice(0, Limit), `${all.length - Limit} more`] : all;
}
