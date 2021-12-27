/*
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import properties from './properties'
import operators from './operators'
import keywords from './keywords'
import functions from './functions'

const docs: string[] = [
    'VMD',
    '============',
    '--------------------------------',
    ''
];

docs.push(`## Properties\n\n`);
docs.push('--------------------------------\n');
for (const name in properties) {
    if (properties[name].isUnsupported) continue

    const names = [name]
    if (properties[name].abbr) names.push(...properties[name].abbr!)
    docs.push(`\`\`\`\n${names.join(', ')}\n\`\`\`\n`);

    if (properties[name]['@desc']) {
        docs.push(`*${properties[name]['@desc']}*\n`);
    }
}

docs.push(`## Operators\n\n`);
docs.push('--------------------------------\n');
operators.forEach(o => {
    if (o.isUnsupported) return

    const names = [o.name]
    if (o.abbr) names.push(...o.abbr!)
    docs.push(`\`\`\`\n${names.join(', ')}\n\`\`\`\n`);

    if (o['@desc']) {
        docs.push(`*${o['@desc']}*\n`);
    }
})

docs.push(`## Keywords\n\n`);
docs.push('--------------------------------\n');
for (const name in keywords) {
    if (!keywords[name].map) continue

    const names = [name]
    if (keywords[name].abbr) names.push(...keywords[name].abbr!)
    docs.push(`\`\`\`\n${names.join(', ')}\n\`\`\`\n`);

    if (keywords[name]['@desc']) {
        docs.push(`*${keywords[name]['@desc']}*\n`);
    }
}

docs.push(`## Functions\n\n`);
docs.push('--------------------------------\n');
for (const name in functions) {
    if (!functions[name].map) continue

    const names = [name]
    docs.push(`\`\`\`\n${names.join(', ')}\n\`\`\`\n`);

    if (functions[name]['@desc']) {
        docs.push(`*${functions[name]['@desc']}*\n`);
    }
}

export default docs.join('\n')