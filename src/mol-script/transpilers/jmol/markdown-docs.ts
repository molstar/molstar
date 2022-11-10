/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 *
 * Adapted from MolQL project
 */

import { properties } from './properties';
import { operators } from './operators';
import { keywords } from './keywords';


const _docs: string[] = [
    'Jmol',
    '============',
    '--------------------------------',
    ''
];

_docs.push(`## Properties\n\n`);
_docs.push('--------------------------------\n');
for (const name in properties) {
    if (properties[name].isUnsupported) continue;

    const names = [name];
    if (properties[name].abbr) names.push(...properties[name].abbr!);
    _docs.push(`\`\`\`\n${names.join(', ')}\n\`\`\`\n`);

    if (properties[name]['@desc']) {
        _docs.push(`*${properties[name]['@desc']}*\n`);
    }
}

_docs.push(`## Operators\n\n`);
_docs.push('--------------------------------\n');
operators.forEach(o => {
    if (o.isUnsupported) return;

    const names = [o.name];
    if (o.abbr) names.push(...o.abbr!);
    _docs.push(`\`\`\`\n${names.join(', ')}\n\`\`\`\n`);

    if (o['@desc']) {
        _docs.push(`*${o['@desc']}*\n`);
    }
});

_docs.push(`## Keywords\n\n`);
_docs.push('--------------------------------\n');
for (const name in keywords) {
    if (!keywords[name].map) continue;

    const names = [name];
    if (keywords[name].abbr) names.push(...keywords[name].abbr!);
    _docs.push(`\`\`\`\n${names.join(', ')}\n\`\`\`\n`);

    if (keywords[name]['@desc']) {
        _docs.push(`*${keywords[name]['@desc']}*\n`);
    }
}

export const docs = _docs.join('\n');
