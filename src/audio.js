// audio.js - the second instrument.
//
// The camera measures rotation by geometry: the sensor reads one row at a
// time, so a fast-spinning object is photographed in slices and comes out
// skewed. The microphone measures the same rotation by sound, because a fan
// or a drill emits a tone at the rate its blades or its motor turn.
//
// Two sensors, two unrelated physical channels, one number. Neither is
// derived from the other and they share no code path beyond the arithmetic
// in dsp.js, which is deliberately the SAME arithmetic for both so that an
// agreement between them is an agreement about the world rather than about
// a particular implementation.
//
// This is also the calibration. The camera cannot know its own row time from
// a datasheet, and the browser will not tell us. But if the microphone says
// the drill is turning at f hertz, and the camera sees that same rotation as
// c cycles per row, then the row time falls out: lineTime = c / f. The
// published methods use a bench LED blinking at a known rate. The device
// already contains a reference.

import { hann, spectrum, prominentPeak } from './dsp.js';

/**
 * Open a microphone stream with every "helpful" processor turned off.
 *
 * Browsers default to voice-call processing: echo cancellation, noise
 * suppression and automatic gain control. All three are designed to remove
 * exactly what we want to keep. Noise suppression in particular treats a
 * steady mechanical tone as noise, which is the literal definition of the
 * signal here.
 */
export async function openMicrophone() {
  if (!navigator.mediaDevices?.getUserMedia) {
    throw new Error('This browser has no getUserMedia. Microphone capture is unavailable.');
  }

  const stream = await navigator.mediaDevices.getUserMedia({
    video: false,
    audio: {
      echoCancellation: false,
      noiseSuppression: false,
      autoGainControl: false,
      channelCount: 1,
    },
  });

  const track = stream.getAudioTracks()[0];
  const settings = track.getSettings ? track.getSettings() : {};

  const ctx = new (window.AudioContext || window.webkitAudioContext)();
  if (ctx.state === 'suspended') await ctx.resume();

  const source = ctx.createMediaStreamSource(stream);
  const analyser = ctx.createAnalyser();
  // A long window, because frequency resolution is what we are buying. At a
  // 48 kHz rate, 16384 samples is about 341 ms and just under 3 Hz per bin.
  analyser.fftSize = 16384;
  analyser.smoothingTimeConstant = 0;   // we do our own averaging, visibly
  source.connect(analyser);

  const warnings = [];
  if (settings.noiseSuppression) warnings.push('Noise suppression is on and the browser refused to disable it. A steady mechanical tone is exactly what it removes.');
  if (settings.autoGainControl) warnings.push('Automatic gain control is on, so amplitude is not trustworthy. Frequency still is.');
  if (settings.echoCancellation) warnings.push('Echo cancellation is on and could not be disabled.');

  return { stream, track, ctx, analyser, settings, warnings, sampleRate: ctx.sampleRate };
}

/**
 * One acoustic measurement, through the same pipeline the camera uses.
 *
 * Deliberately the same window, the same transform and the same peak picker
 * as the image path. If the two instruments agree, that agreement is not an
 * artifact of two different algorithms flattering each other.
 *
 * `band` restricts the search to the mechanically plausible range. A room has
 * speech, traffic and mains hum in it, and none of those are the drill.
 */
export function measureTone(analyser, sampleRate, { minHz = 60, maxHz = 8000 } = {}) {
  const n = analyser.fftSize;
  const buf = new Float32Array(n);
  analyser.getFloatTimeDomainData(buf);

  // Bail out on silence rather than reporting the spectrum of the noise floor.
  let rms = 0;
  for (let i = 0; i < n; i++) rms += buf[i] * buf[i];
  rms = Math.sqrt(rms / n);
  if (rms < 1e-4) return { silent: true, rms };

  const x = new Float64Array(n);
  for (let i = 0; i < n; i++) x[i] = buf[i];
  hann(x);

  const { mag, nfft } = spectrum(x);
  const hzPerBin = sampleRate / nfft;
  const loBin = Math.max(2, Math.floor(minHz / hzPerBin));
  const hiBin = Math.min(mag.length - 2, Math.ceil(maxHz / hzPerBin));

  // Reuse the image path's prominence test, restricted to the audio band.
  const windowed = mag.subarray(0, hiBin + 2);
  const peak = prominentPeak(windowed, nfft, { loBin, baselineHalfWidth: 40 });
  if (!peak) return { silent: false, rms, peak: null };

  return {
    silent: false,
    rms,
    hz: peak.cyclesPerSample * sampleRate,
    prominence: peak.prominence,
    bin: peak.bin,
    hzPerBin,
  };
}

/**
 * Solve the camera's per-row readout time from an acoustic reference.
 *
 * The camera reports a spatial frequency in cycles per row. The microphone
 * reports the same rotation as a temporal frequency in hertz. The row time is
 * the ratio, and it is a property of the sensor rather than of the scene, so
 * once solved it holds for every later measurement on that device.
 *
 *     cyclesPerRow [cycles/row] = f [cycles/s] * lineTime [s/row]
 *
 * Returns seconds per row, or throws rather than returning a number it cannot
 * justify.
 */
export function lineTimeFromAcousticReference(cyclesPerRow, acousticHz) {
  if (!(acousticHz > 0)) throw new Error('lineTimeFromAcousticReference: need a positive acoustic frequency');
  if (!(cyclesPerRow > 0)) throw new Error('lineTimeFromAcousticReference: need a positive cycles/row');
  return cyclesPerRow / acousticHz;
}

/**
 * Do the two instruments agree, and is the implied row time physically real?
 *
 * The second half matters as much as the first. Two wrong numbers can agree.
 * Consumer rolling-shutter sensors read a row in roughly 5 to 60 microseconds,
 * so a solved row time outside that range means the two instruments have
 * locked onto different things (a harmonic, an alias, or the room) even if
 * their ratio looks tidy.
 */
export function crossCheck(cameraCyclesPerRow, acousticHz, { tolPct = 10, minLineTimeUs = 4, maxLineTimeUs = 80 } = {}) {
  if (!(cameraCyclesPerRow > 0) || !(acousticHz > 0)) {
    return { agree: false, reason: 'one instrument has no reading' };
  }
  const lineTime = cameraCyclesPerRow / acousticHz;
  const us = lineTime * 1e6;

  if (us < minLineTimeUs || us > maxLineTimeUs) {
    return {
      agree: false,
      lineTime,
      lineTimeUs: us,
      reason: `the implied row time is ${us.toFixed(1)} us, outside the ${minLineTimeUs} to ${maxLineTimeUs} us range a consumer sensor can physically have. The two instruments are looking at different things, probably a harmonic.`,
    };
  }

  // Try the obvious harmonic relationships. A fan with N blades sounds at N
  // times its shaft rate, and a camera may lock onto either.
  const ratios = [1, 2, 3, 4, 0.5, 1 / 3, 1 / 4];
  const candidates = ratios
    .map((r) => ({ r, us: (cameraCyclesPerRow / (acousticHz * r)) * 1e6 }))
    .filter((c) => c.us >= minLineTimeUs && c.us <= maxLineTimeUs);

  return {
    agree: true,
    lineTime,
    lineTimeUs: us,
    rowRateHz: 1 / lineTime,
    harmonicCandidates: candidates,
    reason: `the implied row time is ${us.toFixed(2)} us, which is physically plausible for a consumer sensor.`,
  };
}
