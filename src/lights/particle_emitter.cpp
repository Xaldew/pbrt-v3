/*
 * Copyright (c) 2019-2023, ARM Limited. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "particle_emitter.h"
#include <cstddef>
#include <vector>
#include "geometry.h"
#include "light.h"
#include "lights/diffuse.h"
#include "pbrt.h"
#include "samplers/random.h"
#include "shapes/triangle.h"
#include "core/paramset.h"
#include "core/sampling.h"
#include "core/reflection.h"
#include "core/scene.h"
#include "spectrum.h"
#include "transform.h"


namespace pbrt
{

const Float Eta_air = 1.000293;


ParticleEmitter::ParticleEmitter(const Transform &LightToWorld,
                                 const MediumInterface &mediumInterface,
                                 int nSamples,
                                 const std::shared_ptr<Shape> &shape,
                                 int nParticles,
                                 Float velocity,
                                 Float range,
                                 Float cherenkov_scale,
                                 Float uniform_scale,
                                 bool randomize,
                                 int seed,
                                 Float ambient_eta)
    : AreaLight(LightToWorld, mediumInterface, nSamples),
      shape(shape),
      cherenkov_scale(cherenkov_scale),
      uniform_scale(uniform_scale),
      area(shape->Area()),
      ambient_eta(ambient_eta)
{
    // Sample the shape and create new particles at each such location.
    RNG rng(seed);
    const Float electron_mass = 9.10938356e-31;
    for (size_t i = 0; i < nParticles; ++i)
    {
        Float pdf;
        Point2f u0 = Point2f{rng.UniformFloat(), rng.UniformFloat()};
        Interaction is = shape->Sample(u0, &pdf);
        Float m = electron_mass;
        float v = velocity;
        Ray r(is.p, Vector3f(is.n), range);
        particles.push_back(Particle{m, v, r});
    }

    if (WorldToLight.HasScale() &&
        dynamic_cast<const Triangle *>(shape.get()) == nullptr)
    {
        Warning(
            "Scaling detected in world to light transformation! "
            "The system has numerous assumptions, implicit and explicit, "
            "that this transform will have no scale factors in it. "
            "Proceed at your own risk; your image may have errors.");
    }
}

ParticleEmitter::~ParticleEmitter() noexcept {}


void ParticleEmitter::Preprocess(const Scene &scene)
{
    MemoryArena arena;
    // Note: All surfaces are assumed to be closed.
    Float sl_length = 0.f;
    Float tot_length = 0.f;
    for (size_t i = 0; i < particles.size(); ++i)
    {
        std::vector<Vector3f> ints;
        std::vector<Spectrum> ft;
        SurfaceInteraction si;
        const auto &p = particles[i];
        RayDifferential r(p.r);
        Float n = ambient_eta;         // Starting refractive_index.
        Float t = 0.;                  // Previous interval start parameter.
        Point3f start(p.r.o);          // Previous interval start point.
        const Primitive *lp = nullptr; // Last primitive hit.
        while (scene.Intersect(r, &si) && p.r.tMax - t > 0.)
        {
            // Is the particle superluminal in the interval [start, si.p]?
            Float nt = (start - si.p).Length();
            if (p.v > 1. / n)
            {
                ints.emplace_back(t, t + nt, n);
                ft.push_back(Frank_Tamm(cherenkov_scale, p.v, n));
                sl_length += nt;
            }
            si.ComputeScatteringFunctions(r, arena);

            // Update interval end-point and IoR of the material we are
            // entering.
            r = si.SpawnRay(p.r.d);
            r.tMax = p.r.tMax - t;
            t += nt;

            // Exiting the primitive, reset IoR to scene IoR.
            if (si.primitive == lp)
            {
                n = ambient_eta;
                lp = nullptr;
            }
            else if (si.bsdf)
            {
                n = si.bsdf->eta;
                lp = si.primitive;
            }
        }
        for (size_t j = 0; j < ints.size(); j++)
        {
            std::cout << "SL interval: " << ints[j] << std::endl;
        }
        ft_spectra.push_back(ft);
        intervals.push_back(ints);
        tot_length += p.r.tMax;
    }
    particle_length = tot_length;
    superluminal_length = sl_length;
    Prob_sl = sl_length / tot_length;
}


Spectrum ParticleEmitter::L(const Interaction &, const Vector3f &) const
{
    // Direct hit on the light source does not actually yield extra light - Only
    // the non-intersectable particles emit light.
    return Spectrum(0.);
}

Spectrum ParticleEmitter::Power() const
{
    // Sum up the emitted light power emitted in each each superluminal
    // interval.
    Spectrum s;
    for (size_t i = 0; i < intervals.size(); ++i)
    {
        for (size_t j = 0; j < intervals[i].size(); ++j)
        {
            auto it = intervals[i][j];
            Float l = it.y - it.x;
            s += ft_spectra[i][j] * l;
        }
    }
    return s;
}


Spectrum ParticleEmitter::Sample_Li(const Interaction &ref,
                                    const Point2f &u,
                                    Vector3f *wi,
                                    Float *pdf,
                                    VisibilityTester *vis) const
{
    // Not applicable for our paper - Thus same as diffuse.cpp.
    ProfilePhase _(Prof::LightSample);
    Interaction pShape = shape->Sample(ref, u, pdf);
    pShape.mediumInterface = mediumInterface;
    if (*pdf == 0 || (pShape.p - ref.p).LengthSquared() == 0) {
        *pdf = 0;
        return 0.f;
    }
    *wi = Normalize(pShape.p - ref.p);
    *vis = VisibilityTester(ref, pShape);
    return L(pShape, -*wi);
}

Float ParticleEmitter::Pdf_Li(const Interaction &ref, const Vector3f &wi) const
{
    // Not applicable for our paper - Thus same as diffuse.cpp.
    ProfilePhase _(Prof::LightPdf);
    return shape->Pdf(ref, wi);
}



Spectrum ParticleEmitter::Sample_Le(const Point2f &u1, const Point2f &u2,
                                    Float time, Ray *ray,
                                    Normal3f *nLight,
                                    Float *pdfPos, Float *pdfDir) const
{
    ProfilePhase _(Prof::LightSample);

    // Sample the light according to our paper.
    if (particles.size() == 0)
    {
        // Sample as a normal area light instead of as a particle emitter.
        Interaction pShape = shape->Sample(u1, pdfPos);
        pShape.mediumInterface = mediumInterface;
        *nLight = pShape.n;

        // Sample a cosine-weighted outgoing direction _w_ for area light.
        Vector3f w = CosineSampleHemisphere(u2);
        *pdfDir = CosineHemispherePdf(w.z);

        Vector3f v1, v2, n(pShape.n);
        CoordinateSystem(n, &v1, &v2);
        w = w.x * v1 + w.y * v2 + w.z * n;
        *ray = pShape.SpawnRay(w);
        return L(pShape, w);
    }
    else
    {
        // u1.x -> Select particle. u1.y -> Select t value for the ray.
        size_t pi = (particles.size() - 1) * u1.x;
        Particle p = particles[pi];
        Float range = p.r.tMax;
        Float t = range * u1.y;

        // Sample a point on the particles path.
        Point3f pp = p.r(t);

        // Is the particle superluminal at the sampled point?
        Spectrum s;
        bool superluminal = false;
        Float n = 1.0;
        for (size_t i = 0; i < intervals[pi].size(); ++i)
        {
            const auto &it = intervals[pi][i];
            if (it.x <= t && t <= it.y)
            {
                superluminal = true;
                n = it.z;
                s = ft_spectra[pi][i];
                break;
            }
        }

        if (superluminal)
        {
            // Select a direction on the Cherenkov emission angle.
            Vector3f w = CherenkovDirection(p.r.d, u2, p.v, n);
            *nLight = (Normal3f)w;
            *ray = Ray(pp, w, Infinity, time, mediumInterface.inside);
        }
        else
        {
            // Select a uniformly spherical direction.
            Vector3f w = UniformSampleSphere(u2);
            *nLight = (Normal3f)w;
            *ray = Ray(pp, w, Infinity, time, mediumInterface.inside);
            s = Spectrum(uniform_scale);
        }

        // Compute the PDFs using a mixture model based on the probability of
        // the sampled point being superluminal.
        *pdfPos = 1. / particle_length;
        *pdfDir = (1. + Prob_sl) * Inv4Pi;

        return s;
    }
}


void ParticleEmitter::Pdf_Le(const Ray &ray, const Normal3f &n,
                             Float *pdfPos, Float *pdfDir) const
{
    ProfilePhase _(Prof::LightPdf);
    Interaction it(ray.o, n, Vector3f(), Vector3f(n), ray.time,
                   mediumInterface);
    *pdfPos = 1 / particle_length;
    *pdfDir = CosineHemispherePdf(Dot(n, ray.d));
}


Vector3f ParticleEmitter::CherenkovDirection(
    const Vector3f &d,
    const Point2f &u,
    Float v,
    Float n) const
{
    Float phi = 2. * Pi * u.x;
    Float theta = std::acos(1. / (v * n));

    const Vector3f ax{1.f, 0.f, 0.f};
    const Vector3f ay{0.f, 1.f, 0.f};

    const Vector3f a = std::fabs(d.x) < 0.1f ? ax : ay;
    const Vector3f r = Normalize(Cross(d, a));

    const Transform R_theta = Rotate(Degrees(theta), r);
    const Transform R_phi = Rotate(Degrees(phi), d);

    return R_phi(R_theta(d));
}


Spectrum ParticleEmitter::Frank_Tamm(Float scale, Float v, Float n) const
{
    // Compute the Frank-Tamm spectrum, and normalize the samples based on the
    // largest (i.e., the first) sample.
    Float lambda[nCIESamples];
    Float vals[nCIESamples];
    const Float a = 0.0072973525693;
    const Float wbeg = CIE_lambda[0] * 1e-9;
    const Float nc = 2 * Pi * a * (1 - (1 / (v * v * n * n))) * (1 / (wbeg * wbeg));
    for (size_t i = 0; i < nCIESamples; ++i)
    {
        Float w = CIE_lambda[i] * 1e-9;
        vals[i] = 2 * Pi * a * (1 - (1 / (v * v * n * n))) * (1 / (w * w));
        vals[i] /= nc;
        vals[i] *= scale;
        lambda[i] = CIE_lambda[i];
    }
    return Spectrum::FromSampled(lambda, vals, nCIESamples);
}


std::shared_ptr<AreaLight>
CreateParticleEmitter(
    const Transform &light2world,
    const Medium *medium,
    const ParamSet &paramSet,
    const std::shared_ptr<Shape> &shape)
{
    int nSamples = paramSet.FindOneInt("samples", paramSet.FindOneInt("nsamples", 1));
    int nParticles = paramSet.FindOneInt("nparticles", 1);
    Float velocity = paramSet.FindOneFloat("velocity", 0.8);
    Float range = paramSet.FindOneFloat("range", 50.);
    Float cherenkov_scale = paramSet.FindOneFloat("cherenkovscale", 1000.f);
    Float uniform_scale = paramSet.FindOneFloat("uniformscale", 50.f);
    bool randomize = paramSet.FindOneBool("randomize", false);
    int seed = paramSet.FindOneInt("seed", 0);
    Float ambient_eta = paramSet.FindOneFloat("ambienteta", Eta_air);
    return std::make_shared<ParticleEmitter>(light2world, medium, nSamples, shape,
                                             nParticles, velocity, range,
                                             cherenkov_scale,
                                             uniform_scale,
                                             randomize,
                                             seed,
                                             ambient_eta);
}


}
