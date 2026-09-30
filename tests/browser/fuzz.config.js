const os = require('node:os');
const path = require('node:path');
const { defineConfig } = require('@playwright/test');
const output = process.env.CABILN_FUZZ_ARTIFACTS || path.join(os.tmpdir(), 'cabiln-browser-fuzz');

module.exports = defineConfig({
  testDir: __dirname,
  testMatch: 'fuzz.spec.cjs',
  timeout: process.env.CABILN_FUZZ_PROFILE === 'deep' ? 900_000 : 180_000,
  expect: { timeout: 10_000 },
  workers: 1,
  retries: 0,
  maxFailures: 1,
  outputDir: path.join(output, 'browser'),
  reporter: [['list'], ['json', { outputFile: path.join(output, 'results.json') }]],
  use: {
    headless: true, channel: process.env.CABILN_BROWSER_CHANNEL || undefined,
    viewport: { width: 1440, height: 1000 }, locale: 'en-GB', colorScheme: 'light',
    // The property owns a fresh context per history and retains its last failed
    // shrink explicitly; runner tracing would double-start those contexts.
    trace: 'off', screenshot: 'only-on-failure',
  },
});
