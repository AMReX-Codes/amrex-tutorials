
#include <MyParticleContainer.H>

namespace amrex {

void 
MyParticleContainer::InitPachinko (std::string initial_particle_file, Real zlen)
{
    BL_PROFILE("MyParticleContainer::InitPachinko()");

    InitFromAsciiFile(initial_particle_file,0);

    int lev = 0;
    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
       int grid = kv.first.first;
       auto& pbox = kv.second.GetArrayOfStructs();
       const int n = pbox.size();

       for (int i = 0; i < n; i++)
       {
            ParticleType& p = pbox[i];

            p.rdata(PIdx::vx) =  0.0;
            p.rdata(PIdx::vy) = -1.0;
#if (AMREX_SPACEDIM == 3)
            p.rdata(PIdx::vz) =  0.0;
#endif

            // We over-write the z-locations to make sure they're in the domain
            p.pos(2) = 0.5 * zlen;
       }
    }
}

void 
MyParticleContainer::AdvectPachinko (Real dt, amrex::Vector<amrex::RealArray>& obstacle_center,
                                     Real obstacle_radius, Real particle_radius)
{
    BL_PROFILE("MyParticleContainer::AdvectPachinko()");

    int lev = 0;

    // Particle radius
    Real prad = particle_radius;

    // Obstracle radiu
    Real crad = obstacle_radius;

    Real rad_squared = (prad+crad)*(prad+crad);

    Real grav = -40.;

    Real restitution_coeff = 0.9;

    const auto prob_lo = Geom(0).ProbLoArray();
    const auto prob_hi = Geom(0).ProbHiArray();

    int num_obstacles = obstacle_center.size();

    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
       int grid = kv.first.first;
       auto& pbox = kv.second.GetArrayOfStructs();
       const int n = pbox.size();

       // std::cout << "Number of particles: " << n << std::endl;

       // NOTE: we assume that all particle motion occurs in a plane!
       // Even if we run in 3-d we assume no motion in the z-direction

       for (int i = 0; i < n; i++)
       {
          ParticleType& p = pbox[i];

          // Let particles stop at the bottom
          if (p.pos(1) >= prob_lo[1] + 0.05) 
          {
            p.pos(0) += 0.5 * dt * p.rdata(PIdx::vx);
            p.pos(1) += 0.5 * dt * p.rdata(PIdx::vy);

            // Accleration under graivty
            p.rdata(PIdx::vy) += dt * grav;

            p.pos(0) += 0.5 * dt * p.rdata(PIdx::vx);
            p.pos(1) += 0.5 * dt * p.rdata(PIdx::vy);

            p.pos(0) += dt * p.rdata(PIdx::vx);
            p.pos(1) += dt * p.rdata(PIdx::vy);

            for (int ind = 0; ind < num_obstacles; ind++)
            {
               Real x_diff = p.pos(0) - obstacle_center[ind][0];
               Real y_diff = p.pos(1) - obstacle_center[ind][1];
               Real diff_sq = x_diff * x_diff + y_diff * y_diff;

               if ( diff_sq < rad_squared )
               {
                 Real diff      = std::sqrt(diff_sq);
                 Real overshoot = (prad+crad) - diff;

                 Real norm_x = x_diff / diff;
                 Real norm_y = y_diff / diff;

                 Real tang_x = -norm_y;
                 Real tang_y =  norm_x;

                 // Incoming velocity dot normal  = (norm_x, norm_y)
                 Real vel_norm = p.rdata(PIdx::vx) * norm_x + 
                                 p.rdata(PIdx::vy) * norm_y;

                 // Incoming velocity dot tangent = (x_tang, y_tang) = (-y_norm, x_norm)
                 Real vel_tang = p.rdata(PIdx::vx) * norm_y + 
                                 p.rdata(PIdx::vy) * norm_x;

                 // Original velocity was (vel_norm) * (norm_x, norm_y)
                 //                    +  (vel_tang) * (tang_x, tang_y)

                 // New      velocity  is MINUS (vel_norm) * (norm_x, norm_y)
                 //                          +  (vel_tang) * (tang_x, tang_y)

                 p.rdata(PIdx::vx) = -vel_norm * norm_x + vel_tang * tang_x;
                 p.rdata(PIdx::vy) = -vel_norm * norm_y + vel_tang * tang_y;

                 p.rdata(PIdx::vx) *= restitution_coeff;
                 p.rdata(PIdx::vy) *= restitution_coeff;

                 // Reflect particle position as well
                 Real ref_pos_x = obstacle_center[ind][0] + (prad+crad)*norm_x;
                 Real ref_pos_y = obstacle_center[ind][1] + (prad+crad)*norm_y;

                 p.pos(0) = ref_pos_x + overshoot * norm_x;
                 p.pos(1) = ref_pos_y + overshoot * norm_y;
               }
            }

            // Bounce off left wall
            if (p.pos(0) < (prob_lo[0]+prad)) 
            {
               p.pos(0) = 2.0*(prob_lo[0]+prad) - p.pos(0);
               p.rdata(PIdx::vx) = -p.rdata(PIdx::vx);
               p.rdata(PIdx::vx) *= restitution_coeff;
               p.rdata(PIdx::vy) *= restitution_coeff;
            }

            // Bounce off right wall
            if (p.pos(0) > (prob_hi[0]-prad)) 
            {
               p.pos(0) = 2.0*(prob_hi[0]-prad) - p.pos(0);
               p.rdata(PIdx::vx) = -0.9*p.rdata(PIdx::vx);
               p.rdata(PIdx::vx) *= restitution_coeff;
               p.rdata(PIdx::vy) *= restitution_coeff;
            }
         }
       }
    }
}
}
