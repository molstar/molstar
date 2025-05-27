const { createApp, createExample } = require('./webpack.config.common.js');

const examples = [
    'proteopedia-wrapper',
    'basic-wrapper',
    'lighting',
    'alpha-orbitals',
    'mvs-stories',
    'ihm-restraints',
    'interactions',
    'ligand-editor'
];

module.exports = [
    createApp('viewer', 'molstar'),
    createApp('docking-viewer', 'molstar'),
    createApp('mesoscale-explorer', 'molstar'),
    createApp('mvs-stories', 'mvsStories', { filename: 'mvs-stories.js', cssFilename: 'mvs-stories.css' }),
    ...examples.map(createExample)
];