import { defineConfig } from "eslint/config";
import globals from "globals";
import typescriptEslint from "@typescript-eslint/eslint-plugin";
import tsParser from "@typescript-eslint/parser";

export default defineConfig([{
    ignores: [
        "node_modules/*",
        "build/*",
        "deploy/*",
        "docs/site/*",
        "lib/*",
        "eslint.config.mjs",
        "build.mjs",
    ]
},{
    languageOptions: {
        globals: {
            ...globals.browser,
            ...globals.node,
        },

        ecmaVersion: 2018,
        sourceType: "module",

        parserOptions: {
            ecmaFeatures: {
                impliedStrict: true,
            },
        },
    },

    rules: {
        indent: "off",
        "arrow-parens": ["off", "as-needed"],
        "brace-style": ["error", "1tbs", {
            allowSingleLine: true,
        }],
        "comma-spacing": "off",
        "space-infix-ops": "off",
        "comma-dangle": "off",
        eqeqeq: ["error", "smart"],
        "import/order": "off",
        "no-eval": "warn",
        "no-extend-native": "warn",
        "no-new-wrappers": "warn",
        "no-trailing-spaces": "error",
        "no-unsafe-finally": "warn",
        "no-self-compare": "warn",
        "no-var": "error",
        "spaced-comment": "error",
        semi: "warn",
        "no-restricted-syntax": ["error", {
            selector: "ExportDefaultDeclaration",
            message: "Default exports are not allowed",
        }],
        "no-throw-literal": "error",
        "key-spacing": "error",
        "object-curly-spacing": ["error", "always"],
        "array-bracket-spacing": "error",
        "space-in-parens": "error",
        "computed-property-spacing": "error",
        "prefer-const": ["error", {
            destructuring: "all",
            ignoreReadBeforeAssign: false,
        }],
        "space-before-function-paren": "off",
        "func-call-spacing": "off",
        "no-multi-spaces": "error",
        "block-spacing": "error",
        "keyword-spacing": "warn",
        "space-before-blocks": "error",
        "semi-spacing": "error",
        "no-constant-binary-expression": "error",
    },
}, {
    files: ["**/*.ts", "**/*.tsx"],

    plugins: {
        "@typescript-eslint": typescriptEslint,
    },

    languageOptions: {
        parser: tsParser,
        ecmaVersion: 5,
        sourceType: "module",

        parserOptions: {
            project: ["tsconfig.eslint.json"],
        },
    },

    rules: {
        "@typescript-eslint/ban-types": "off",
        "@typescript-eslint/class-name-casing": "off",        
        "@typescript-eslint/member-delimiter-style": ["off", {
            multiline: {
                delimiter: "none",
                requireLast: true,
            },

            singleline: {
                delimiter: "semi",
                requireLast: false,
            },
        }],
        "@typescript-eslint/prefer-namespace-keyword": "warn",
        "@typescript-eslint/semi": ["off", null],
    },
}]);