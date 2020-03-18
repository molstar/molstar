const common = require('./webpack.config.common.js');
const createApp = common.createApp; 
const createEntry = common.createEntry;
const createBrowserTest = common.createBrowserTest;
const createNodeApp = common.createNodeApp;

module.exports = [
    createApp('viewer'),
    createApp('basic-wrapper'),
    createEntry('examples/proteopedia-wrapper/index', 'examples/proteopedia-wrapper', 'index'),
    createEntry('apps/demos/lighting/index', 'demos/lighting', 'index'),
    createNodeApp('state-docs'),
    createApp('model-server-query'),

    createBrowserTest('font-atlas'),
    createBrowserTest('marching-cubes'),
    createBrowserTest('render-lines'),
    createBrowserTest('render-mesh'),
    createBrowserTest('render-shape'),
    createBrowserTest('render-spheres'),
    createBrowserTest('render-structure'),
    createBrowserTest('render-text'),
]