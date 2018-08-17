const { generate } = require('graphql-code-generator')
const path = require('path')

const basePath = path.join(__dirname, '..', '..', 'src', 'mol-model-props', 'rcsb', 'graphql')

generate({
    args: [
        path.join(basePath, 'symmetry.gql.ts')
    ],
    schema: 'http://rest-experimental.rcsb.org/graphql',
    template: 'typescript',
    out: path.join(basePath),
    skipSchema: true,
    overwrite: true,
    config: path.join(__dirname, 'codegen.json')
}, true).then(
    console.log('done')
).catch(e => console.error(e))