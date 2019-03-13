const { parse } = require('graphql');
const { readFileSync } = require('fs');

module.exports = function(docString, config) {
    const str = readFileSync(docString, { encoding: 'utf-8' }).trim()
                    .replace(/^export default `/, '')
                    .replace(/`$/, '')
    return [
        {
            filePath: docString,
            content: parse(str)
        }
    ];
};