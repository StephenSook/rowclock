// Run with:  node --test tests/
//
// These tests check the solver against closed-form answers we can write down
// in advance. The synthetic signals here are a UNIT TEST of an algorithm, not
// evidence about the product: nothing in this file is ever quoted as a
// measurement of a real scene. Every product claim is measured on a real
// camera frame and is reported separately.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  fft,
  nextPow2,
  detrend,
  spectrum,
  dominantPeak,
  analyseProfile,
  lineTimeFromKnownSource,
  hzFromCyclesPerRow,
  aliasNote,
} from '../src/dsp.js';

/** A sinusoid at a known cycles-per-sample rate, optionally over a scene ramp. */
function synth(n, cyclesPerSample, { amp = 1, phase = 0, dc = 0, ramp = 0, noise = 0, seed = 1 } = {}) {
  // Deterministic LCG so a failure is reproducible.
  let s = seed >>> 0;
  const rnd = () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296) - 0.5;
  const x = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    x[i] = dc
      + ramp * (i / n)
      + amp * Math.sin(2 * Math.PI * cyclesPerSample * i + phase)
      + noise * rnd();
  }
  return x;
}

test('nextPow2 rounds up to a power of two', () => {
  assert.equal(nextPow2(1), 1);
  assert.equal(nextPow2(1080), 2048);
  assert.equal(nextPow2(2048), 2048);
  assert.equal(nextPow2(2049), 4096);
});

test('fft rejects a non-power-of-two length rather than returning nonsense', () => {
  assert.throws(() => fft(new Float64Array(6), new Float64Array(6)), /power of two/);
});

test('fft of a unit impulse is flat at magnitude 1', () => {
  const n = 64;
  const re = new Float64Array(n);
  const im = new Float64Array(n);
  re[0] = 1;
  fft(re, im);
  for (let k = 0; k < n; k++) {
    assert.ok(Math.abs(Math.hypot(re[k], im[k]) - 1) < 1e-12, `bin ${k} was not 1`);
  }
});

test('fft puts a pure tone in exactly the bin it belongs in', () => {
  const n = 256;
  const k0 = 17;
  const re = new Float64Array(n);
  const im = new Float64Array(n);
  for (let i = 0; i < n; i++) re[i] = Math.cos((2 * Math.PI * k0 * i) / n);
  fft(re, im);
  let best = 1;
  for (let k = 1; k < n / 2; k++) if (Math.hypot(re[k], im[k]) > Math.hypot(re[best], im[best])) best = k;
  assert.equal(best, k0);
});

test('detrend removes DC and a scene ramp but keeps the tone', () => {
  const n = 1080;
  const cps = 0.05;
  const clean = synth(n, cps);
  const dirty = synth(n, cps, { dc: 120, ramp: 60 });

  const dt = detrend(dirty, 65);

  // A centered moving-average high-pass cannot drive the mean to exactly zero,
  // because the window is truncated at the two edges of the array and is
  // therefore asymmetric there. For a ramp of height `ramp` the edge residual
  // is of order ramp*half/n, so asserting 1e-6 asserts a perfection this
  // filter structurally cannot deliver. The property that actually matters is
  // that the residual is negligible NEXT TO THE TONE we are trying to keep.
  const mean = dt.reduce((a, b) => a + b, 0) / n;
  const toneAmplitude = 1; // synth default
  assert.ok(
    Math.abs(mean) < 0.01 * toneAmplitude,
    `residual DC ${mean} is not small relative to the tone`,
  );

  // The property a high-pass actually guarantees is spectral, not statistical:
  // the scene's energy lives in the lowest bins, and after detrending those
  // bins must be small next to the tone's own bin. Comparing means is the
  // wrong instrument here, because the interior legitimately still contains
  // the tone we set out to keep, so its mean is not a measure of cleanliness.
  const before = spectrum(Float64Array.from(dirty));
  const after = spectrum(dt);
  const toneBin = Math.round(cps * after.nfft);

  const lowBandEnergy = (s) => {
    let e = 0;
    for (let k = 1; k < 4; k++) e += s.mag[k] * s.mag[k];
    return e;
  };

  const ratioBefore = lowBandEnergy(before) / (before.mag[toneBin] ** 2);
  const ratioAfter = lowBandEnergy(after) / (after.mag[toneBin] ** 2);
  assert.ok(
    ratioAfter < ratioBefore / 100,
    `detrend only cut low-band energy from ${ratioBefore.toExponential(2)} to ${ratioAfter.toExponential(2)}`,
  );

  // The tone should survive with most of its amplitude intact.
  const pClean = analyseProfile(clean);
  const pDirty = analyseProfile(dirty);
  assert.ok(pClean && pDirty);
  assert.ok(
    Math.abs(pClean.cyclesPerSample - pDirty.cyclesPerSample) < 1e-4,
    `scene content moved the estimate: ${pClean.cyclesPerSample} vs ${pDirty.cyclesPerSample}`,
  );
});

test('dominantPeak recovers a known frequency to better than a tenth of a bin', () => {
  // Deliberately off-bin so the parabolic interpolation is doing real work.
  const n = 1080;
  const nfft = nextPow2(n); // 2048
  const trueCps = 137.4 / nfft;
  const x = synth(n, trueCps, { dc: 100, ramp: 40, noise: 0.02 });

  const peak = analyseProfile(x);
  assert.ok(peak, 'no peak found');

  const errBins = Math.abs(peak.bin - 137.4);
  assert.ok(errBins < 0.1, `peak off by ${errBins.toFixed(3)} bins (got ${peak.bin.toFixed(3)})`);
  assert.ok(peak.snr > 5, `snr too low: ${peak.snr.toFixed(1)}`);
});

test('interpolation beats raw bin picking on an off-bin tone', () => {
  const n = 1024;
  const nfft = nextPow2(n);
  const trueBin = 80.5; // exactly between two bins, the worst case
  const x = synth(n, trueBin / nfft);
  const { mag, nfft: nf } = spectrum(x);
  const peak = dominantPeak(mag, nf, 4);

  let rawBest = 4;
  for (let k = 4; k < mag.length - 1; k++) if (mag[k] > mag[rawBest]) rawBest = k;

  const rawErr = Math.abs(rawBest - trueBin);
  const interpErr = Math.abs(peak.bin - trueBin);
  assert.ok(interpErr <= rawErr, `interpolation made it worse: ${interpErr} vs ${rawErr}`);
});

test('a flat profile with no tone does not manufacture a confident peak', () => {
  const n = 1080;
  const x = synth(n, 0, { amp: 0, dc: 128, noise: 0.5, seed: 7 });
  const peak = analyseProfile(x);
  // A peak object may still come back; what must not happen is a high SNR.
  if (peak) assert.ok(peak.snr < 5, `manufactured a peak at snr ${peak.snr.toFixed(1)}`);
});

test('calibration round-trips: known source in, correct unknown out', () => {
  // Pretend the sensor reads a row every 20 microseconds (50 kHz row rate).
  const trueLineTime = 20e-6;

  // A 120 Hz mains lamp therefore bands at:
  const calCyclesPerRow = 120 * trueLineTime; // 0.0024 cycles/row

  const solved = lineTimeFromKnownSource(calCyclesPerRow, 120);
  assert.ok(Math.abs(solved - trueLineTime) < 1e-12, `line time off: ${solved}`);

  // Now an unknown 400 Hz source through the same sensor.
  const unknownCyclesPerRow = 400 * trueLineTime;
  const hz = hzFromCyclesPerRow(unknownCyclesPerRow, solved);
  assert.ok(Math.abs(hz - 400) < 1e-9, `recovered ${hz} Hz, expected 400`);
});

test('calibration refuses impossible inputs instead of returning a number', () => {
  assert.throws(() => lineTimeFromKnownSource(0.002, 0), /positive/);
  assert.throws(() => lineTimeFromKnownSource(0, 120), /positive/);
  assert.throws(() => hzFromCyclesPerRow(0.002, 0), /positive/);
});

test('end to end on a synthetic 400 Hz scene, recovered within 1 percent', () => {
  const trueLineTime = 20e-6;
  const trueHz = 400;
  const rows = 1080;

  // Build the row profile a camera would actually produce: scene content,
  // vignetting, sensor noise, and the illumination modulation on top.
  const cps = trueHz * trueLineTime;
  const profile = synth(rows, cps, { amp: 3.5, dc: 110, ramp: -35, noise: 1.2, seed: 42 });

  const peak = analyseProfile(profile);
  assert.ok(peak, 'no peak');
  const hz = hzFromCyclesPerRow(peak.cyclesPerSample, trueLineTime);

  const pctErr = (Math.abs(hz - trueHz) / trueHz) * 100;
  assert.ok(pctErr < 1, `recovered ${hz.toFixed(1)} Hz, ${pctErr.toFixed(2)} percent error`);
});

test('aliasNote reports the fold point and the ambiguous twin', () => {
  const lineTime = 20e-6;        // 50 kHz row rate
  const note = aliasNote(400, lineTime);
  assert.ok(Math.abs(note.rowRateHz - 50000) < 1e-6);
  assert.ok(Math.abs(note.nyquistHz - 25000) < 1e-6);
  assert.equal(note.unambiguous, true);
  assert.ok(note.aliases.some((f) => Math.abs(f - 49600) < 1e-6), 'missing the lower alias');

  const folded = aliasNote(30000, lineTime);
  assert.equal(folded.unambiguous, false, 'a frequency above Nyquist must be flagged');
});
