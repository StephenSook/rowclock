// dsp.js - the numerical core. No dependencies, no build step, no network.
//
// Everything here operates on plain Float64Array. It is written to be callable
// identically from the browser page and from a Node harness, so the same code
// produces the numbers in the live demo and the numbers in the tests.

/**
 * In-place iterative radix-2 Cooley-Tukey FFT.
 * re and im must be the same length and that length must be a power of two.
 * Decimation in time, with the bit-reversal permutation done up front.
 */
export function fft(re, im) {
  const n = re.length;
  if (n !== im.length) throw new Error(`fft: re/im length mismatch (${n} vs ${im.length})`);
  if (n < 2 || (n & (n - 1)) !== 0) throw new Error(`fft: length must be a power of two, got ${n}`);

  // Bit-reversal permutation.
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      let t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }

  // Butterflies.
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * Math.PI) / len;
    const wRe = Math.cos(ang);
    const wIm = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let curRe = 1, curIm = 0;
      for (let k = 0; k < len / 2; k++) {
        const aRe = re[i + k],            aIm = im[i + k];
        const bRe = re[i + k + len / 2],  bIm = im[i + k + len / 2];
        const tRe = bRe * curRe - bIm * curIm;
        const tIm = bRe * curIm + bIm * curRe;
        re[i + k] = aRe + tRe;  im[i + k] = aIm + tIm;
        re[i + k + len / 2] = aRe - tRe;
        im[i + k + len / 2] = aIm - tIm;
        const nextRe = curRe * wRe - curIm * wIm;
        curIm = curRe * wIm + curIm * wRe;
        curRe = nextRe;
      }
    }
  }
}

/** Next power of two >= n. */
export function nextPow2(n) {
  let p = 1;
  while (p < n) p <<= 1;
  return p;
}

/**
 * Remove DC and low-order scene structure.
 *
 * This matters more than the FFT does. A row-mean profile of a real photograph
 * is dominated by the SCENE (a bright ceiling, a dark desk, lens vignetting),
 * and that content sits at very low spatial frequency. Subtracting a centered
 * moving average of width `win` is a crude high-pass that leaves the flicker
 * band untouched while removing the picture.
 *
 * Returns a new Float64Array; does not modify the input.
 */
export function detrend(x, win = 65) {
  const n = x.length;
  const out = new Float64Array(n);
  if (win < 3) { out.set(x); return out; }
  const half = win >> 1;

  // Prefix sums so the moving average is O(n) rather than O(n*win).
  const pre = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) pre[i + 1] = pre[i] + x[i];

  for (let i = 0; i < n; i++) {
    const a = Math.max(0, i - half);
    const b = Math.min(n, i + half + 1);
    const mean = (pre[b] - pre[a]) / (b - a);
    out[i] = x[i] - mean;
  }
  return out;
}

/** Hann window, applied in place. Reduces spectral leakage from the finite record. */
export function hann(x) {
  const n = x.length;
  for (let i = 0; i < n; i++) {
    x[i] *= 0.5 * (1 - Math.cos((2 * Math.PI * i) / (n - 1)));
  }
  return x;
}

/**
 * Magnitude spectrum of a real signal, zero-padded to the next power of two.
 * Returns { mag, nfft } where mag has nfft/2 usable bins.
 * Bin k corresponds to k / nfft cycles per sample.
 */
export function spectrum(x) {
  const nfft = nextPow2(x.length);
  const re = new Float64Array(nfft);
  const im = new Float64Array(nfft);
  re.set(x);
  fft(re, im);
  const half = nfft >> 1;
  const mag = new Float64Array(half);
  for (let k = 0; k < half; k++) {
    mag[k] = Math.hypot(re[k], im[k]);
  }
  return { mag, nfft };
}

/**
 * Locate the dominant peak and refine it to sub-bin precision.
 *
 * Quadratic interpolation over the log-magnitudes of the peak bin and its two
 * neighbours. Without this the frequency resolution is one bin, which at
 * nfft=2048 over 1080 rows is a coarse answer; with it the estimate is good to
 * a small fraction of a bin, which is what makes a three-significant-figure
 * readout honest rather than decorative.
 *
 * `loBin` skips the near-DC bins that survive detrending.
 * Returns { bin, cyclesPerSample, magnitude, snr } or null if nothing stands out.
 */
export function dominantPeak(mag, nfft, loBin = 4) {
  const n = mag.length;
  if (loBin >= n - 1) return null;

  let best = loBin;
  for (let k = loBin; k < n - 1; k++) {
    if (mag[k] > mag[best]) best = k;
  }
  if (best <= 0 || best >= n - 1) return null;
  if (!(mag[best] > 0)) return null;

  // Quadratic (parabolic) interpolation in the log domain.
  const eps = 1e-12;
  const a = Math.log(mag[best - 1] + eps);
  const b = Math.log(mag[best] + eps);
  const c = Math.log(mag[best + 1] + eps);
  const denom = a - 2 * b + c;
  const delta = denom === 0 ? 0 : (0.5 * (a - c)) / denom;
  const refined = best + Math.max(-0.5, Math.min(0.5, delta));

  // Median of the whole band as a noise floor, so the ratio is robust to the
  // peak itself and to a handful of other strong bins.
  const sorted = Float64Array.from(mag.subarray(loBin)).sort();
  const median = sorted[sorted.length >> 1] || eps;

  return {
    bin: refined,
    cyclesPerSample: refined / nfft,
    magnitude: mag[best],
    snr: mag[best] / median,
  };
}

/**
 * The whole pipeline for one profile: detrend, window, transform, pick the peak.
 * `profile` is a per-row (or per-sample) series. Returns the peak, or null.
 */
export function analyseProfile(profile, { detrendWin = 65, loBin = 4 } = {}) {
  if (!profile || profile.length < 32) return null;
  const d = hann(detrend(profile, detrendWin));
  const { mag, nfft } = spectrum(d);
  return dominantPeak(mag, nfft, loBin);
}

/**
 * Solve the sensor's per-row readout time from a source of KNOWN frequency.
 *
 * This is the calibration that makes the whole instrument possible, and it is
 * why no datasheet is needed. A mains-powered lamp modulates at twice the line
 * frequency (120 Hz on a 60 Hz supply) because instantaneous power goes as
 * v(t)^2 and therefore peaks twice per cycle. Photograph one, measure the
 * banding period in CYCLES PER ROW, and the row period falls out:
 *
 *     cyclesPerRow [cycles/row] = f_known [cycles/s] * lineTime [s/row]
 *  => lineTime = cyclesPerRow / f_known
 *
 * Returns seconds per row.
 */
export function lineTimeFromKnownSource(cyclesPerRow, knownHz) {
  if (!(knownHz > 0)) throw new Error('lineTimeFromKnownSource: knownHz must be positive');
  if (!(cyclesPerRow > 0)) throw new Error('lineTimeFromKnownSource: cyclesPerRow must be positive');
  return cyclesPerRow / knownHz;
}

/** Convert a measured spatial frequency to a temporal one, given the calibration. */
export function hzFromCyclesPerRow(cyclesPerRow, lineTimeSeconds) {
  if (!(lineTimeSeconds > 0)) throw new Error('hzFromCyclesPerRow: lineTime must be positive');
  return cyclesPerRow / lineTimeSeconds;
}

/**
 * The frequency the rolling shutter CANNOT distinguish from f.
 *
 * Row sampling is still sampling, so it aliases. With a row rate of
 * 1/lineTime, anything above half that folds back. Reporting the fold points
 * alongside a reading is the difference between an instrument and a number
 * generator, and it is the first thing an imaging engineer will ask about.
 */
export function aliasNote(hz, lineTimeSeconds) {
  const rowRate = 1 / lineTimeSeconds;
  const nyquist = rowRate / 2;
  return {
    rowRateHz: rowRate,
    nyquistHz: nyquist,
    unambiguous: hz < nyquist,
    // The other candidates that produce an identical banding period.
    aliases: [rowRate - hz, rowRate + hz].filter((f) => f > 0 && f < rowRate * 3),
  };
}
