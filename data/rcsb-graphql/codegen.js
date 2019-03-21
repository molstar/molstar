const { generate } = require('graphql-code-generator')
const path = require('path')

const basePath = path.join(__dirname, '..', '..', 'src', 'mol-model-props', 'rcsb', 'graphql')

generate({
    schema: 'http://rest-dev.rcsb.org/graphql',
    documents: {
        [path.join(basePath, 'symmetry.gql.ts')]: {
            loader: path.join(__dirname, 'loader.js')
        },
    },
    generates: {
        [path.join(basePath, 'types.ts')]: {
            plugins: ['time', 'typescript-common', 'typescript-client']
        }
    },
    overwrite: true,
    config: path.join(__dirname, 'codegen.json')
}, true).then(
    () => console.log('done')
).catch(e => console.error(e))