const { createApp, createExample } = require('./webpack.config.common.js');

const examples = ['proteopedia-wrapper', 'basic-wrapper', 'lighting'];

module.exports = [
    createApp('viewer', 'molstar'),
    createApp('docking-viewer', 'molstar'),
    ...examples.map(createExample)
]