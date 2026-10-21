#!/usr/bin/env node
/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import '@molstar/plugin/state/transforms/catalog';
import { StateTransformer, StateObject } from '@molstar/core/state';
import { StringBuilder } from '@molstar/core/util';
import * as fs from 'fs';
import { paramsToMd } from '@molstar/state-docs-cli/pd-to-md';
import type { PluginSpec } from '@molstar/plugin/spec';

const withParticles = (spec: PluginSpec): PluginSpec => ({
  ...spec,
  registry: [...(spec.registry ?? []), ...ParticlesRegistry],
});
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { ParamDefinition } from '@molstar/core/util/param-definition';

// The default plugin and the opt-in particles, so the docs list every built-in representation and theme
const ctx = new PluginContext(withParticles(DefaultPluginSpec()));
await ctx.init();

const builder = StringBuilder.create();

function typeToString(o: StateObject.Ctor[]) {
  if (o.length === 0) return '()';
  return o.map((o) => o.name).join(' | ');
}

function writeTransformer(t: StateTransformer) {
  StringBuilder.write(
    builder,
    `## <a name="${t.id.replace('.', '-')}"></a>${t.id} :: ${typeToString(t.definition.from)} -> ${typeToString(t.definition.to)}`,
  );
  StringBuilder.newline(builder);
  if (t.definition.display.description) {
    StringBuilder.write(builder, `*${t.definition.display.description}*`);
    StringBuilder.newline(builder);
  }
  StringBuilder.newline(builder);
  if (t.definition.params) {
    const params = t.definition.params(void 0, ctx);
    StringBuilder.write(builder, `### Parameters`);
    StringBuilder.newline(builder);
    StringBuilder.write(builder, paramsToMd(params));
    StringBuilder.newline(builder);

    StringBuilder.write(builder, `### Default Parameters`);
    StringBuilder.newline(builder);
    StringBuilder.write(
      builder,
      `\`\`\`js\n${JSON.stringify(ParamDefinition.getDefaultValues(params), null, 2)}\n\`\`\``,
    );
    StringBuilder.newline(builder);
  }
  StringBuilder.write(builder, '----------------------------');
  StringBuilder.newline(builder);
}

const transformers = StateTransformer.getAll();

StringBuilder.write(builder, '# Mol* Plugin State Transformer Reference');
StringBuilder.newline(builder);
StringBuilder.newline(builder);
transformers.forEach((t) => {
  StringBuilder.write(builder, `* [${t.id}](#${t.id.replace('.', '-')})`);
  StringBuilder.newline(builder);
});
StringBuilder.newline(builder);
StringBuilder.write(builder, '----------------------------');
StringBuilder.newline(builder);
transformers.forEach((t) => writeTransformer(t));

fs.mkdirSync('docs/state', { recursive: true });
fs.writeFileSync(`docs/state/transforms.md`, StringBuilder.getString(builder));
