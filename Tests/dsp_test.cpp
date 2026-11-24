// Standalone checks for the Angel Shift DSP core (no JUCE needed).
//   g++ -O2 -std=c++17 Tests/dsp_test.cpp -o dsp_test && ./dsp_test
#include "../Source/FormantShifter.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>
#include <complex>

static const double kPi = 3.14159265358979323846;
static const double FS  = 48000.0;
static int failures = 0;

#define CHECK(cond, ...) do { if (!(cond)) { ++failures; std::printf("  FAIL: "); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// ---- helpers -------------------------------------------------------------
static std::vector<float> makeVowel (double seconds, double f0, const double F[3], const double BW[3])
{
    const int N = (int) (seconds * FS);
    std::vector<double> x (N, 0.0);
    const double period = FS / f0;
    for (double p = 0; p < N; p += period) x[(int) p] = 1.0;
    // glottal tilt
    double s = 0; for (auto& v : x) { s = 0.9 * s + v; v = s; }
    for (int f = 0; f < 3; ++f)
    {
        const double r = std::exp (-kPi * BW[f] / FS), th = 2 * kPi * F[f] / FS;
        const double a1 = -2 * r * std::cos (th), a2 = r * r;
        double y1 = 0, y2 = 0;
        for (auto& v : x) { const double y = v - a1 * y1 - a2 * y2; y2 = y1; y1 = y; v = y; }
    }
    double pk = 0; for (double v : x) pk = std::max (pk, std::abs (v));
    std::vector<float> out (N);
    for (int i = 0; i < N; ++i) out[i] = (float) (0.5 * x[i] / pk);
    return out;
}

static double rms (const std::vector<float>& x, int a, int b)
{
    double s = 0; for (int i = a; i < b; ++i) s += (double) x[i] * x[i];
    return std::sqrt (s / std::max (1, b - a));
}

static double estimateF0 (const std::vector<float>& x, int a, int b)
{
    const int lmin = (int) (FS / 400), lmax = (int) (FS / 60);
    std::vector<double> c (lmax + 1, 0.0); double best = -1e30;
    for (int l = lmin; l <= lmax; ++l)
    {
        double s = 0; for (int i = a; i + l < b; ++i) s += (double) x[i] * x[i + l];
        c[l] = s; best = std::max (best, s);
    }
    for (int l = lmin + 1; l < lmax; ++l)            // first strong local maximum (avoids picking 2T)
        if (c[l] >= 0.9 * best && c[l] >= c[l - 1] && c[l] >= c[l + 1]) return FS / l;
    return 0.0;
}

// independent formant estimator: LPC order 12 + Durand-Kerner roots
[[maybe_unused]] static void estimateFormants (const std::vector<float>& x, int a, int b, double out[3])
{
    const int p = 12, N = b - a;
    std::vector<double> w (N);
    double prev = 0;
    for (int i = 0; i < N; ++i)
    {
        const double v = x[a + i];
        w[i] = (v - 0.95 * prev) * (0.5 - 0.5 * std::cos (2 * kPi * i / (N - 1)));
        prev = v;
    }
    double r[p + 1];
    for (int k = 0; k <= p; ++k) { r[k] = 0; for (int i = k; i < N; ++i) r[k] += w[i] * w[i - k]; }
    r[0] *= 1.0001;
    double A[p + 1] = { 1 }, err = r[0];
    for (int i = 1; i <= p; ++i)
    {
        double acc = r[i]; for (int j = 1; j < i; ++j) acc += A[j] * r[i - j];
        const double k = -acc / err; double t[p + 1]; std::copy (A, A + p + 1, t);
        for (int j = 1; j < i; ++j) A[j] = t[j] + k * t[i - j];
        A[i] = k; err *= 1 - k * k;
    }
    using C = std::complex<double>;
    C z[p];
    for (int i = 0; i < p; ++i) z[i] = std::polar (0.9, 2 * kPi * i / p + 0.3);
    for (int it = 0; it < 500; ++it)
        for (int i = 0; i < p; ++i)
        {
            C pv (1, 0), den (1, 0);
            for (int k = 1; k <= p; ++k) pv = pv * z[i] + A[k];
            for (int j = 0; j < p; ++j) if (j != i) den *= (z[i] - z[j]);
            z[i] -= pv / den;
        }
    std::vector<double> f;
    for (int i = 0; i < p; ++i)
        if (z[i].imag() > 0)
        {
            const double fr = std::arg (z[i]) * FS / (2 * kPi), bw = -std::log (std::abs (z[i])) * FS / kPi;
            if (fr > 200 && fr < 6000 && bw < 500) f.push_back (fr);
        }
    std::sort (f.begin(), f.end());
    for (int i = 0; i < 3; ++i) out[i] = i < (int) f.size() ? f[i] : 0.0;
}


static std::vector<double> harmonicEnvelope (const std::vector<float>& x, int a, int b, double f0)
{
    std::vector<double> E;
    for (int k = 1; k * f0 < 7000; ++k)
    {
        double re = 0, im = 0; const double w = 2 * kPi * k * f0 / FS;
        for (int i = a; i < b; ++i) { re += x[i] * std::cos (w * i); im -= x[i] * std::sin (w * i); }
        E.push_back (20 * std::log10 (std::sqrt (re * re + im * im) / (b - a) + 1e-9));
    }
    return E;
}
static double envAt (const std::vector<double>& E, double f0, double f)
{
    const double pos = f / f0 - 1.0;
    if (pos <= 0) return E.front();
    const int i = (int) pos; if (i + 1 >= (int) E.size()) return E.back();
    const double t = pos - i; return E[i] * (1 - t) + E[i + 1] * t;
}


// ---- tiny FFT + PSD for the formant test ------------------------------------
static void fft (std::vector<std::complex<double>>& a)
{
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) { size_t bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) std::swap (a[i], a[j]); }
    for (size_t len = 2; len <= n; len <<= 1)
    {
        const double ang = -2 * kPi / len; const std::complex<double> wl (std::cos (ang), std::sin (ang));
        for (size_t i = 0; i < n; i += len) { std::complex<double> w (1); for (size_t j = 0; j < len / 2; ++j) { auto u = a[i + j], v = a[i + j + len / 2] * w; a[i + j] = u + v; a[i + j + len / 2] = u - v; w *= wl; } }
    }
}
static std::vector<double> psd (const std::vector<float>& x, int a, int b)   // 4096-pt Welch, Hann
{
    const int N = 4096; std::vector<double> P (N / 2, 0.0); int cnt = 0;
    for (int s = a; s + N <= b; s += N / 2)
    {
        std::vector<std::complex<double>> f (N);
        for (int i = 0; i < N; ++i) f[i] = x[s + i] * (0.5 - 0.5 * std::cos (2 * kPi * i / N));
        fft (f); for (int i = 0; i < N / 2; ++i) P[i] += std::norm (f[i]); ++cnt;
    }
    for (auto& v : P) v /= cnt; return P;
}
static double peakIn (const std::vector<double>& P, double f0, double f1)   // smoothed peak frequency in [f0,f1]
{
    const double df = FS / 4096; double best = -1; double bf = f0;
    for (int i = (int) (f0 / df) + 3; i < (int) (f1 / df) - 3; ++i)
    { double s = 0; for (int k = -3; k <= 3; ++k) s += P[i + k]; if (s > best) { best = s; bf = i * df; } }
    return bf;
}
static std::vector<float> noiseVowel (double seconds, const double F[3], const double BW[3])
{
    std::mt19937 g (11); std::normal_distribution<double> nd (0, 1);
    const int N = (int) (seconds * FS); std::vector<double> x (N);
    double s = 0; for (auto& v : x) { s = 0.9 * s + nd (g); v = s; }
    for (int f = 0; f < 3; ++f)
    {
        const double r = std::exp (-kPi * BW[f] / FS), th = 2 * kPi * F[f] / FS, a1 = -2 * r * std::cos (th), a2 = r * r;
        double y1 = 0, y2 = 0; for (auto& v : x) { const double y = v - a1 * y1 - a2 * y2; y2 = y1; y1 = y; v = y; }
    }
    double pk = 0; for (double v : x) pk = std::max (pk, std::abs (v));
    std::vector<float> o (N); for (int i = 0; i < N; ++i) o[i] = (float) (0.4 * x[i] / pk); return o;
}

static std::vector<float> run (const std::vector<float>& in, float st, int block = 512)
{
    angel::FormantShifter fx;
    fx.setShiftSemitones (st);
    fx.prepare (FS);
    std::vector<float> L = in, R = in;
    for (size_t i = 0; i < L.size(); i += block)
        fx.process (L.data() + i, R.data() + i, (int) std::min<size_t> (block, L.size() - i));
    return L;
}

int main()
{
    const double F[3]  = { 700, 1200, 2600 }, BW[3] = { 90, 100, 140 };
    const double F2[3] = { 330, 2000, 2900 }, BW2[3] = { 70, 110, 150 };   // an "ee"
    auto vowel = makeVowel (3.0, 100.0, F, BW);

    // 1. identity at 0 st --------------------------------------------------
    {
        std::printf ("[1] identity at 0 st\n");
        std::mt19937 g (1); std::normal_distribution<float> nd (0.f, 0.05f);
        auto x = vowel; for (auto& v : x) v += nd (g);
        auto y = run (x, 0.f);
        double md = 0; for (size_t i = 0; i < x.size(); ++i) md = std::max (md, (double) std::abs (x[i] - y[i]));
        std::printf ("  max |y-x| = %.3g\n", md);
        CHECK (md < 1e-4, "not an identity (%g)", md);
    }

    // 2. formants move, pitch stays, level stays --------------------------
    for (int vi = 0; vi < 2; ++vi)
    {
        auto x = vi == 0 ? vowel : makeVowel (3.0, 200.0, F2, BW2);
        const double f0in = vi == 0 ? 100.0 : 200.0;
        const int per = (int) std::lround (FS / f0in);
        const int a = per * (int) (1.5 * FS / per), b = a + per * 40;
        std::printf ("[2] vowel %d (f0 = %.0f Hz)\n", vi, f0in);
        for (float st : { -12.f, -6.f, -3.f, 3.f, 6.f, 12.f })
        {
            auto y = run (x, st);
            const double f0 = estimateF0 (y, a, a + 8192);
            const double dB = 20 * std::log10 (rms (y, a, b) / rms (x, a, b));
            std::printf ("    %+5.0f st: f0 = %.1f Hz  level %+.1f dB\n", st, f0, dB);
            CHECK (std::abs (f0 / f0in - 1.0) < 0.01, "pitch changed: %.2f", f0);
            CHECK (std::abs (dB) < 4.0, "level jump %.1f dB", dB);
        }
    }

    // 2b. formant positions measured on a noise-excited vowel ------------------
    {
        std::printf ("[2b] formant peaks (noise-excited vowel, F1/F2 = 700/1500 Hz)\n");
        const double FN[3] = { 700, 1500, 2800 }, BN[3] = { 110, 130, 180 };
        auto x = noiseVowel (6.0, FN, BN);
        const int a = (int) (2.0 * FS), b = (int) (6.0 * FS);
        auto P0 = psd (x, a, b);
        std::printf ("    input: F1 %.0f  F2 %.0f\n", peakIn (P0, 450, 1100), peakIn (P0, 1150, 2200));
        for (float st : { -6.f, -3.f, 3.f, 6.f })
        {
            auto y = run (x, st);
            auto P = psd (y, a, b);
            const double ratio = std::pow (2.0, st / 12.0);
            const double f1 = peakIn (P, 700 * ratio * 0.7, 700 * ratio * 1.35), f2 = peakIn (P, 1500 * ratio * 0.8, 1500 * ratio * 1.25);
            std::printf ("    %+3.0f st: F1 %.0f (target %.0f)   F2 %.0f (target %.0f)\n", st, f1, 700 * ratio, f2, 1500 * ratio);
            CHECK (std::abs (f1 / (700 * ratio) - 1) < 0.15, "F1 off target");
            CHECK (std::abs (f2 / (1500 * ratio) - 1) < 0.15, "F2 off target");
        }
    }

    // 3. stability fuzz ----------------------------------------------------
    {
        std::printf ("[3] stability fuzz\n");
        std::mt19937 g (7); std::uniform_real_distribution<float> ur (-1.f, 1.f);
        std::uniform_int_distribution<int> blk (1, 700);
        std::uniform_real_distribution<float> sh (-12.f, 12.f);
        const int N = (int) (12 * FS);
        std::vector<float> x (N);
        for (int i = 0; i < N; ++i)
        {
            const int seg = (int) (i / FS) % 6;
            float v = 0;
            switch (seg)
            {
                case 0: v = ur (g); break;                                   // full-scale noise
                case 1: v = vowel[i % vowel.size()]; break;                  // vowel
                case 2: v = ((i / 60) & 1) ? 1.f : -1.f; break;              // square (clipped)
                case 3: v = 0.f; break;                                      // silence
                case 4: v = (i % 4800 == 0) ? 1.f : 0.f; break;              // impulses
                case 5: v = 0.9f; break;                                     // DC
            }
            x[i] = v;
        }
        angel::FormantShifter fx; fx.prepare (FS);
        std::vector<float> L = x, R = x;
        double pk = 0; bool finite = true;
        for (int i = 0; i < N;)
        {
            const int n = std::min (blk (g), N - i);
            fx.setShiftSemitones (sh (g));
            fx.process (L.data() + i, R.data() + i, n);
            i += n;
        }
        for (int i = 0; i < N; ++i) { if (!std::isfinite (L[i]) || !std::isfinite (R[i])) finite = false; pk = std::max (pk, (double) std::abs (L[i])); }
        std::printf ("  finite = %d, peak = %.2f\n", (int) finite, pk);
        CHECK (finite, "NaN/Inf in output");
        CHECK (pk < 8.0, "output peak too large: %.2f", pk);
    }

    // 4. transients: timing + level (0.5 ms RMS envelope) --------------------
    {
        std::printf ("[4] transients\n");
        std::mt19937 g (3); std::normal_distribution<float> nd (0.f, 1.f);
        for (float st : { -6.f, 6.f })
        {
            auto x = vowel;
            const int burst = (int) (1.2 * FS), click = (int) (2.0 * FS);
            for (int i = 0; i < (int) (0.01 * FS); ++i) x[burst + i] += 0.8f * nd (g) * std::exp (-i / (0.002f * (float) FS));
            x[click] += 0.8f;
            auto y = run (x, st);
            auto env = [&] (const std::vector<float>& v, int a, int b, int& at)
            {
                const int W = 24; double best = 0; at = a;
                for (int i = a; i < b; ++i) { double s = 0; for (int k = 0; k < W; ++k) s += (double) v[i + k] * v[i + k]; s = std::sqrt (s / W); if (s > best) { best = s; at = i; } }
                return best;
            };
            int ai, ao;
            const double ei = env (x, click - 100, click + 100, ai), eo = env (y, click - 100, click + 100, ao);
            int bi, bo;
            const double fi = env (x, burst, burst + 480, bi), fo = env (y, burst, burst + 480, bo);
            std::printf ("  %+3.0f st: click env %+.1f dB (moved %+d smp)   burst env %+.1f dB (moved %+d smp)\n",
                         st, 20 * std::log10 (eo / ei), ao - ai, 20 * std::log10 (fo / fi), bo - bi);
            CHECK (std::abs (ao - ai) <= 48 && std::abs (bo - bi) <= 48, "transient moved by more than 1 ms");
            CHECK (std::abs (20 * std::log10 (eo / ei)) < 6.0, "click level changed");
            CHECK (std::abs (20 * std::log10 (fo / fi)) < 5.0, "burst level changed");
        }
    }

    // 5. sibilants / noise are left alone ----------------------------------
    {
        std::printf ("[5] sibilant protection\n");
        std::mt19937 g (5); std::normal_distribution<float> nd (0.f, 0.2f);
        std::vector<float> x ((int) (2 * FS)); float p = 0;
        for (auto& v : x) { const float n = nd (g); v = n - p; p = n; }            // high-passed noise ("s")
        for (float st : { -8.f, 8.f })
        {
            auto y = run (x, st);
            const int a = (int) (1 * FS), b = (int) (2 * FS);
            double e = 0, s = 0; for (int i = a; i < b; ++i) { e += std::pow ((double) y[i] - x[i], 2); s += (double) x[i] * x[i]; }
            const double dB = 10 * std::log10 (e / s + 1e-20);
            std::printf ("  %+3.0f st: error vs. dry = %.1f dB\n", st, dB);
            CHECK (dB < -12.0, "noise got shifted (%.1f dB)", dB);
        }
    }

    // 6. CPU ----------------------------------------------------------------
    {
        angel::FormantShifter fx; fx.prepare (FS); fx.setShiftSemitones (5.f);
        std::vector<float> L = vowel, R = vowel; L.resize ((size_t) (FS * 3)); R = L;
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i + 512 <= L.size(); i += 512) fx.process (L.data() + i, R.data() + i, 512);
        const double sec = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
        std::printf ("[6] CPU: %.2f%% of one core (stereo, 48 kHz, this machine)\n", 100.0 * sec / 3.0);
    }

    std::printf (failures ? "\n%d CHECK(S) FAILED\n" : "\nALL CHECKS PASSED\n", failures);
    return failures ? 1 : 0;
}
