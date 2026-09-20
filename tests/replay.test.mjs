// The replay layer turns a recording into a pointing outcome. These tests
// check the arithmetic against answers written down in advance, and pin the
// honesty properties: the generous metric, the swept parameter, and the
// clinical grading.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {
  WCAG_MIN_PX, WCAG_ENHANCED_PX,
  cpiToPxPerMm, recordingToPath, holdFraction, targetNeverLeft,
  targetForHoldRate, hitCurve, sampleAt, updrsBand,
} from '../src/replay.js';

const REC = new URL('../src/data/006_StretchHold_RightWrist.txt', import.meta.url);

test('the WCAG constants are the published ones', () => {
  assert.equal(WCAG_MIN_PX, 24);        // SC 2.5.8, Level AA
  assert.equal(WCAG_ENHANCED_PX, 44);   // SC 2.5.5, Level AAA
});

test('cpi to px/mm is exact and refuses nonsense', () => {
  assert.ok(Math.abs(cpiToPxPerMm(800) - 800 / 25.4) < 1e-12);
  assert.ok(Math.abs(cpiToPxPerMm(25.4) - 1) < 1e-12);
  assert.throws(() => cpiToPxPerMm(0), /positive/);
});

test('a synthetic circular path gives the hold fraction geometry predicts', () => {
  // Closed form for a cursor tracing a circle of radius R against a square
  // target of half-side h, centred together. The cursor is inside when
  // |cos a| <= h/R AND |sin a| <= h/R. Squaring and adding those gives
  // 1 <= 2(h/R)^2, so NO point of the circle is inside unless h >= R/sqrt(2).
  //
  // The first version of this test asserted "about a third" for h = R/2,
  // which was arithmetic done in the head and never checked. The true answer
  // is exactly zero, and the code was right.
  const n = 3600, R = 10;
  const x = new Float64Array(n), y = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const a = (2 * Math.PI * i) / n;
    x[i] = R * Math.cos(a); y[i] = R * Math.sin(a);
  }
  const path = { x, y, n, fs: 100, seconds: n / 100 };

  // Side 2R, half-side R: every point qualifies.
  assert.ok(Math.abs(holdFraction(path, 1, 2 * R + 0.01) - 1) < 1e-9, 'a big enough square must hold all of it');

  // Side R, half-side R/2, which is below R/sqrt(2) = 0.707R: exactly none.
  assert.equal(holdFraction(path, 1, R), 0);

  // Just below the threshold: still none. Just above: some, but not all.
  const hThresh = R / Math.SQRT2;
  assert.equal(holdFraction(path, 1, 2 * hThresh * 0.99), 0);
  const partial = holdFraction(path, 1, 2 * 0.9 * R);
  assert.ok(partial > 0 && partial < 1, `expected a partial hold, got ${partial}`);

  assert.ok(Math.abs(targetNeverLeft(path, 1) - 2 * R) < 1e-6);
});

test('holdFraction refuses a non-positive target', () => {
  const path = { x: new Float64Array([0]), y: new Float64Array([0]), n: 1, fs: 100, seconds: 0.01 };
  assert.throws(() => holdFraction(path, 1, 0), /positive/);
});

test('targetForHoldRate brackets the size that achieves a given hold rate', () => {
  const n = 2000, R = 8;
  const x = new Float64Array(n), y = new Float64Array(n);
  for (let i = 0; i < n; i++) { const a = (2 * Math.PI * i) / n; x[i] = R * Math.cos(a); y[i] = R * Math.sin(a); }
  const path = { x, y, n, fs: 100, seconds: n / 100 };

  const s = targetForHoldRate(path, 1, 0.95);
  assert.ok(s, 'should find a size');
  assert.ok(holdFraction(path, 1, s) >= 0.95, 'the returned size must actually achieve the rate');
  assert.ok(holdFraction(path, 1, s * 0.9) < 0.95, 'and it must be near the boundary, not oversized');
});

test('REAL DATA: the committed Parkinson\'s recording reproduces the headline', () => {
  const path = recordingToPath(fs.readFileSync(REC, 'utf8'));

  // Frequency must land in the published Parkinsonian rest-tremor band,
  // 4 to 7 Hz (Bhatia et al. 2018, Movement Disorders consensus).
  assert.ok(path.hz > 4 && path.hz < 7, `tremor frequency ${path.hz?.toFixed(2)} Hz is outside the published PD band`);
  assert.ok(path.prominence > 10, `prominence ${path.prominence?.toFixed(1)} is too low to call this a clean tremor`);

  // Amplitude must be clinically plausible: above physiological tremor
  // (0.0973 mm, Duval 2005) and within the MDS-UPDRS ladder.
  assert.ok(path.p2pMm > 0.5, `${path.p2pMm.toFixed(3)} mm is barely above physiological tremor`);
  assert.ok(path.p2pMm < 100, `${path.p2pMm.toFixed(1)} mm is off the clinical scale entirely`);

  // The finding itself, at the common default sensitivity.
  const hold = holdFraction(path, cpiToPxPerMm(800), WCAG_MIN_PX);
  assert.ok(hold < 0.5, `hold rate at the WCAG minimum was ${(hold * 100).toFixed(0)}%, which would refute the headline`);
});

test('REAL DATA: the finding survives the whole sensitivity sweep above 400 cpi', () => {
  const path = recordingToPath(fs.readFileSync(REC, 'utf8'));
  const curve = hitCurve(path);
  for (const row of curve.grid) {
    if (row.cpi < 400) continue;
    const wcagCell = row.cells.find((c) => c.sizePx === WCAG_MIN_PX);
    assert.ok(
      wcagCell.hold < 0.6,
      `at ${row.cpi} cpi the WCAG minimum held ${(wcagCell.hold * 100).toFixed(0)}% of the time, breaking the claim`,
    );
    assert.ok(row.neverLeftPx > WCAG_MIN_PX, `at ${row.cpi} cpi the tremor fits inside a 24 px target`);
  }
});

test('REAL DATA: a healthy control on the same task is not comparable', () => {
  const pd = recordingToPath(fs.readFileSync(REC, 'utf8'));
  const hc = recordingToPath(fs.readFileSync(new URL('../src/data/001_StretchHold_RightWrist.txt', import.meta.url), 'utf8'));
  assert.ok(
    pd.prominence > hc.prominence,
    `the healthy control scored a cleaner tremor line (${hc.prominence?.toFixed(1)}) than the patient (${pd.prominence?.toFixed(1)})`,
  );
});

test('the hold rate is monotonic in target size, as it must be', () => {
  const path = recordingToPath(fs.readFileSync(REC, 'utf8'));
  const ppm = cpiToPxPerMm(800);
  let prev = -1;
  for (const s of [8, 16, 24, 44, 64, 128, 256, 512]) {
    const h = holdFraction(path, ppm, s);
    assert.ok(h >= prev - 1e-12, `hold rate fell from ${prev} to ${h} as the target grew`);
    prev = h;
  }
});

test('sampleAt loops and interpolates rather than stepping', () => {
  const path = recordingToPath(fs.readFileSync(REC, 'utf8'));
  const a = sampleAt(path, 0.0, 1);
  const b = sampleAt(path, path.seconds, 1);   // exactly one loop later
  assert.ok(Math.abs(a.x - b.x) < 1e-9 && Math.abs(a.y - b.y) < 1e-9, 'the loop is not seamless');

  const mid = sampleAt(path, 0.005, 1);        // half a sample in
  assert.ok(Number.isFinite(mid.x) && Number.isFinite(mid.y));
});

test('updrsBand grades honestly against the published ladder', () => {
  assert.equal(updrsBand(3.775).grade, 1);          // our measurement
  assert.equal(updrsBand(3.775).label, 'slight');
  assert.equal(updrsBand(20).grade, 2);             // 2 cm
  assert.equal(updrsBand(50).grade, 3);             // 5 cm
  assert.equal(updrsBand(150).grade, 4);            // 15 cm
});

test('the committed data manifest describes what is actually there', () => {
  const m = JSON.parse(fs.readFileSync(new URL('../src/data/manifest.json', import.meta.url), 'utf8'));
  assert.equal(m.license, 'CC BY-NC-SA 4.0');
  assert.ok(m.url.includes('physionet.org'));
  assert.ok(m.records.length >= 4);
  for (const r of m.records) {
    const p = new URL(`../src/data/${r.file}`, import.meta.url);
    assert.ok(fs.existsSync(p), `manifest lists ${r.file} but it is not committed`);
    assert.ok(r.condition, `${r.file} has no condition label`);
    assert.equal(r.samplingRateHz, 100);
  }
});
