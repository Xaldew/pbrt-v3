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

#ifndef SRC_LIGHTS_PARTICLE_EMITTER_H
#define SRC_LIGHTS_PARTICLE_EMITTER_H

#include "geometry.h"
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#include "pbrt.h"
#include "light.h"
#include "primitive.h"

namespace pbrt {

class ParticleEmitter : public AreaLight
{
public:

    ParticleEmitter(const Transform &LightToWorld,
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
                    Float ambient_ior);

    virtual ~ParticleEmitter() noexcept;

    virtual void Preprocess(const Scene &scene);

    Spectrum L(const Interaction &intr, const Vector3f &w) const;

    Spectrum Power() const;
    Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wo,
                       Float *pdf, VisibilityTester *vis) const;
    Float Pdf_Li(const Interaction &, const Vector3f &) const;
    Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
                       Ray *ray, Normal3f *nLight, Float *pdfPos,
                       Float *pdfDir) const;
    void Pdf_Le(const Ray &,
                const Normal3f &,
                Float *pdfPos,
                Float *pdfDir) const;

protected:

    struct Particle
    {
        Float m; /**< Mass. */
        Float v; /**< Velocity as a percentage of the speed of light. */
        Ray r;   /**< Origin and propagation direction. */
    };

    /** @brief Compute the direction of a Cherenkov photon. */
    Vector3f CherenkovDirection(const Vector3f &d, const Point2f &u, Float v, Float n) const;

    /** @brief Compute the Frank-Tamm spectrum for a particle. */
    std::vector<Float> Frank_Tamm(Float v, Float n) const;

protected:
    std::shared_ptr<Shape> shape;
    std::vector<Particle> particles;
    std::vector<std::vector<Vector3f>> intervals; /**< Superluminal locations. */
    std::vector<std::vector<Spectrum>> ft_spectra; /**< Spectra for the interval. */
    const Float cherenkov_scale;
    const Float uniform_scale;
    const Float ambient_eta;
    const Float area;
    Float particle_length; /**< Total particle lengths. */
    Float superluminal_length; /**< Total superluminal length. */
    Float Prob_sl; /**< Probability of the particle being superluminal. */
};


std::shared_ptr<AreaLight> CreateParticleEmitter(
    const Transform &light2world,
    const Medium *medium,
    const ParamSet &paramSet,
    const std::shared_ptr<Shape> &shape);

}  // namespace pbrt


#endif /* SRC_LIGHTS_PARTICLE_EMITTER_H */
