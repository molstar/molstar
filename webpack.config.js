const { createApp, createExample, createBrowserTest } = require('./webpack.config.common.js');

const apps = ['viewer'];
const examples = ['proteopedia-wrapper', 'basic-wrapper', 'lighting'];
const tests = [
    'font-atlas',
    'marching-cubes',
    'render-lines', 'render-mesh', 'render-shape', 'render-spheres', 'render-structure', 'render-text',
    'parse-xtc'
];

module.exports = [
    ...apps.map(createApp),
    ...examples.map(createExample),
    ...tests.map(createBrowserTest)
]