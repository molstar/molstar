const { createApp, createExample } = require('./webpack.config.common.js');

const apps = ['viewer'];
const examples = ['proteopedia-wrapper', 'basic-wrapper', 'lighting'];

module.exports = [
    ...apps.map(createApp),
    ...examples.map(createExample)
]