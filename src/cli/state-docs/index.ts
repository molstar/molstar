#!/usr/bin/env node
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as _ from '../../mol-plugin-state/transforms';
import { StateTransformer, StateObject } from '../../mol-state';
import { StringBuilder } from '../../mol-util';
import * as fs from 'fs';
import { paramsToMd } from './pd-to-md';
import { PluginContext } from '../../mol-plugin/context';
import { ParamDefinition } from '../../mol-util/param-definition';

// force the transform to be evaluated
_.StateTransforms.Data.Download.id;

// Empty plugin context
const ctx = new PluginContext({
    actions: [],
    behaviors: []
});

const builder = StringBuilder.create();

function typeToString(o: StateObject.Ctor[]) {
    if (o.length === 0) return '()';
    return o.map(o => o.name).join(' | ');
}

function writeTransformer(t: StateTransformer) {
    StringBuilder.write(builder, `## <a name="${t.id.replace('.', '-')}"></a>${t.id} :: ${typeToString(t.definition.from)} -> ${typeToString(t.definition.to)}`);
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
        StringBuilder.write(builder, `\`\`\`js\n${JSON.stringify(ParamDefinition.getDefaultValues(params), null, 2)}\n\`\`\``);
        StringBuilder.newline(builder);
    }
    StringBuilder.write(builder, '----------------------------');
    StringBuilder.newline(builder);
}

const transformers = StateTransformer.getAll();

StringBuilder.write(builder, '# Mol* Plugin State Transformer Reference');
StringBuilder.newline(builder);
StringBuilder.newline(builder);
transformers.forEach(t => {
    StringBuilder.write(builder, `* [${t.id}](#${t.id.replace('.', '-')})`);
    StringBuilder.newline(builder);
});
StringBuilder.newline(builder);
StringBuilder.write(builder, '----------------------------');
StringBuilder.newline(builder);
transformers.forEach(t => writeTransformer(t));

fs.writeFileSync(`docs/state/transforms.md`, StringBuilder.getString(builder));