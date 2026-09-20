# Prior art, and what is actually ours

We did not invent this physics. Using a rolling shutter as a sub-frame temporal
sampler is a published technique with a patent and a decade of literature behind
it, and several of the specific steps in this project are the published steps.
This page exists so that a reader who knows that literature does not have to
catch us claiming otherwise.

Everything below was checked against sources we read. Where a paper was
paywalled and only the abstract was available, it says so.

---

## The method is not ours

**US Patent 11,776,238 B2, US Naval Research Laboratory.** This is the general
method, and it is the closest prior art to what we built. Row-mean intensity as
the observable, an explicit term for the inter-frame blanking interval, and the
Lomb-Scargle periodogram to handle the resulting uneven sampling. Verbatim from
the patent:

> "we can calculate the row time based on the pixel row, image width, and
> horizontal blanking, and can express the time of the sample taken from row m
> of frame n for video with row time r and frame rate v as: n v + m · r"

and

> "using 30 fps video, we correctly identify signals at over 1000 Hz, as well as
> identify signals from complex waveforms including triangle, square, and
> combinations of sinusoids"

**Zebra-stripe angle for shaft speed.** "Rotating machinery speed extraction
through smartphone video acquisition from a radial viewpoint", Mechanical
Systems and Signal Processing, 2023. ABSTRACT ONLY, paywalled. Verbatim:

> "A method is proposed to measure the speed from the deformation of a zebra
> pattern as a consequence of the sequential readout of a rolling shutter
> camera."

and, on beating the frame-rate limit:

> "This characteristic period is the basis for the developed methodology and
> allows measuring constant and varying speeds above the Nyquist limit related
> to the camera's frame rate."

**Smartphone video tachometry.** André, Leclère, Anastasio, Benaïcha, Billon,
Birem, Bonnardot et al., "Using a smartphone camera to analyse rotating and
vibrating systems: Feedback on the SURVISHNO 2019 contest", MSSP 154 (2021)
107553. Paywalled; the open HAL copy was gated. Verbatim from the abstract:

> "the consideration of the unavoidable slight phase shift between the
> acquisition of each pixel opens up the possibility to perform a dynamic
> analysis at frequencies that are much higher than the video frame rate."

**Calibrating the row time by counting bands under a known-frequency light.**
`horshack-dpreview/RollingShutter`. Verbatim:

> "if we know the light source is 60Hz (120 transitions/second) and the sensor
> captures 6 bands, we can calculate the sensor readout time via 1000/120*6"

Our `lineTimeFromKnownSource()` is that same relation stated per row. That
project measures dedicated cameras with a bench Arduino LED at 500 Hz and
publishes no smartphone values; its maintainer has said he has no plans to test
phones.

**Measuring a camera's internal timing from images alone.** Sensors 26(16):5231.
Verbatim:

> "Two internal timing parameters are estimated directly from the captured
> images without access to the internal camera timing: the row readout period
> (34.38 μs), obtained from the spatial periodicity of the stripes, and the
> effective integration time (490 μs)."

**Other rolling-shutter vibrometry we are not claiming to have originated:**
Hong et al., Optics Letters 48(15):3837 (2023), laser speckle, 14.285 kHz
recovered at 70 fps with a 35 µs inter-row period. Sheinin, Chan, O'Toole,
Narasimhan, "Dual-Shutter Optical Vibration Sensing", CVPR 2022, up to 63 kHz
from sensors rated at 130 Hz. Zhao et al., electronic rolling shutter, errors
under 2% from 300 Hz to 20 kHz and under 5% from 40 to 300 Hz.

**Aliasing disambiguation** is the stroboscope canon: start high, step down, and
halve to confirm. Not new.

**Browser rolling-shutter analysis is already commoditized for light flicker.**
FlickerHz and FlickerTest both read row brightness in a web page and classify
against IEEE 1789-2015. We checked this before building and it is why this
project does not do flicker measurement as its product.

---

## What is actually ours

Three things, and they are narrow on purpose.

**1. The implementation is un-shipped.** An adversarial search for a shipped
rotation or vibration tool built on the rolling-shutter sampling mechanism, on
any platform, found none. Every camera tachometer app we found is frame-rate
stroboscopic and therefore capped near the frame rate, and treats rolling
shutter as a defect. The iOS "Video Tachometer" app says so in its own
description:

> "when using a format with better quality, the Rolling Shutter effect is
> stronger. Therefore, the best quality format can be used at speeds less than
> about 2000 rpm."

That is the opposite of using the effect as the sampler. Packaging the academic
method as something a stranger runs in a browser tab with no install and no
account is engineering, not new physics, and we say so.

**2. Calibration from the device's own microphone.** The published methods
calibrate the row time against an external reference at a known frequency: a
bench Arduino LED, or a separate video of a blinking LED. A rotating machine
emits sound at the rate it turns, so the same rotation can be measured
acoustically and optically at once, and the row time is the ratio of the two.
The reference is already inside the laptop. We have not found this in the
literature, and if it is there we will cite it.

**3. A measured account of what the browser camera path does to the signal.**
This is the part the literature does not have. The research pass concluded that
a device-by-device characterization of whether `getUserMedia` preserves the
rolling-shutter signature is undocumented, and our own build turned into exactly
that study by accident. See `docs/FALSE-GREENS.md` for the four convincing
spurious readings this instrument produced on real hardware before it learned to
refuse them.

---

## Constraints we cannot engineer around, and did not try to hide

**`powerLineFrequency` is not a web constraint.** It exists as a 2015 proposal to
the W3C media capture spec that was never adopted, and as an OS-level control
(V4L2 `power_line_frequency`, DirectShow `put_PowerlineFrequency`). Our capture
code requests it because it costs nothing on a platform that might one day
support it, but it cannot work today and is labelled as such in the source.

**Exposure control is unavailable in practice.** `exposureMode` and
`exposureTime` are defined only in the Image Capture spec, and only as optional
advanced constraints. iOS Safari does not support ImageCapture at all. The
MacBook this was developed on reports no exposure control whatsoever, and the
page says so on screen rather than implying the capture path is clean.

**That has a direct physical consequence.** Rolling-shutter banding from a light
only exists if a row's exposure is shorter than one cycle of the source, 8.3 ms
at 120 Hz. An indoor webcam on auto-exposure typically sits at 16 to 33 ms, so
the modulation is integrated to zero before any code sees it. This is why this
project measures rotation, whose observable is geometric, rather than
illumination flicker, whose observable is photometric and which the browser can
erase before we get it.

**Android anti-banding defaults to AUTO** and actively cancels the signal, with
no standardized way to disable it from the web.

---

## Sources

| What | Where | Access |
|---|---|---|
| Row-sampling method, blanking, Lomb-Scargle | US Patent 11,776,238 B2 | full text read |
| Zebra-stripe shaft speed | MSSP, S0888327023007446 | abstract only, paywalled |
| Smartphone video tachometry | MSSP 154 (2021) 107553 | abstract only, paywalled |
| Line-laser rolling-shutter vibrometry, 36 kHz | ResearchGate 403584444 | abstract only |
| Laser speckle, 14.285 kHz at 70 fps | Optics Letters 48(15):3837 | abstract read |
| Dual-shutter, 63 kHz | CVPR 2022, par.nsf.gov/servlets/purl/10395505 | full text read |
| Row-time calibration by band counting | github.com/horshack-dpreview/RollingShutter | full text read |
| Internal timing from images, 34.38 µs | Sensors 26(16):5231 | full text read |
| Browser flicker tools | flickertest.com, github.com/tonym128/flash_light | read |
| Media capture constraints | w3.org/TR/image-capture, w3.org/TR/mediacapture-streams | read |
| Anti-banding defaults to AUTO | Android CameraCharacteristics | read |
