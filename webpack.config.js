const { createApp, createExample, createBrowserTest } = require('./webpack.config.common.js');

const examples = [
    'proteopedia-wrapper',
    'basic-wrapper',
    'lighting',
    'alpha-orbitals',
    'alphafolddb-pae',
    'mvs-stories',
    'ihm-restraints',
    'interactions',
    'ligand-editor'
];
const tests = [
    'font-atlas',
    'marching-cubes',
    'render-lines', 'render-mesh', 'render-shape', 'render-spheres', 'render-structure', 'render-structure-grid', 'render-text',
    'parse-xtc'
];

module.exports = [
    createApp('viewer', 'molstar'),
    createApp('docking-viewer', 'molstar'),
    createApp('mesoscale-explorer', 'molstar'),
    createApp('mvs-stories', 'mvsStories', { filename: 'mvs-stories.js', cssFilename: 'mvs-stories.css' }),
    ...examples.map(createExample),
    ...tests.map(createBrowserTest)
];
