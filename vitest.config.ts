// vitest.config.ts
import { defineConfig } from 'vitest/config';

export const vitestConfig = defineConfig({
    test: {
        globals: true,
        environment: 'jsdom',
    },
});
