/** @type {import('jest').Config} */
module.exports = {
    moduleFileExtensions: ['ts', 'js'],
    transform: {
      '\\.ts$': 'ts-jest'
    },
    moduleDirectories: ['node_modules', 'lib'],
    testEnvironmentOptions: {
      url: 'http://localhost/'
    },
    testRegex: '/__tests__/jest/.*\\.test\\.ts$'
};