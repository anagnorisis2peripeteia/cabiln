const os = require('node:os');
const path = require('node:path');
const { defineConfig } = require('@playwright/test');

// Install with npm ci and npx playwright install chromium in this directory.
// The Python environment must have the project's [web] dependencies installed.
// CABILN_APP_ROOT selects a frozen checkout; CABILN_PYTHON selects its interpreter.
// CABILN_BROWSER_CHANNEL=chrome uses an installed Chrome instead of Chromium.
const outputDir = process.env.CABILN_BROWSER_ARTIFACTS ||
  path.join(os.tmpdir(), 'cabiln-browser-acceptance');

module.exports = defineConfig({
  testDir: __dirname,
  testMatch: '*.spec.js',
  timeout: 60_000,
  expect: { timeout: 15_000 },
  workers: 1,
  retries: 0,
  outputDir,
  reporter: [
    ['list'],
    ['json', { outputFile: path.join(outputDir, 'results.json') }],
  ],
  use: {
    headless: true,
    channel: process.env.CABILN_BROWSER_CHANNEL || undefined,
    viewport: { width: 1440, height: 1000 },
    locale: 'en-GB',
    colorScheme: 'light',
    screenshot: 'only-on-failure',
    trace: 'retain-on-failure',
  },
});
